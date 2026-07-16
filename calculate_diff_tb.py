import os
import sys

# =======================================================
# 🚀 救命参数：破解并发死锁与线程爆炸 (必须放在所有的 import 之前)
# =======================================================
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"  # 禁用 NetCDF 底层文件锁 (解决卡死的核心)
os.environ["OMP_NUM_THREADS"] = "1"            # 限制 OpenMP 线程
os.environ["MKL_NUM_THREADS"] = "1"            # 限制 MKL 线程
os.environ["OPENBLAS_NUM_THREADS"] = "1"       # 限制 OpenBLAS 线程
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"     # 限制 Mac/Linux 向量计算线程
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import glob
import time
import math
import numpy as np
import netCDF4 as nc
import torch
import multiprocessing as mp
import concurrent.futures
from tqdm import tqdm
from soil_rtm import SoilRTM_TB_Calculator

def worker_process(file_idx, file_path, total_files, output_dir):
    """
    单文件处理的独立 worker，自动分配 GPU 并重定向日志
    直接将结果写入独立的 NetCDF 文件。
    """
    # =======================================================
    # 🚀 在子进程内部强行限制 PyTorch 的 CPU 线程，把算力全部让给 GPU
    # =======================================================
    torch.set_num_threads(1)
    
    # 获取进程 ID 用于绑定 GPU 和日志文件
    current_process = mp.current_process()
    worker_id = current_process._identity[0] if current_process._identity else 1
    
    # 日志输出重定向
    log_dir = "log_files_compare"
    os.makedirs(log_dir, exist_ok=True)
    log_file_path = os.path.join(log_dir, f"worker_{worker_id}.log")
    
    log_file = open(log_file_path, 'a', encoding='utf-8', buffering=1)
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = log_file
    sys.stderr = log_file

    try:
        # 动态分配 GPU (利用身份 ID 对 4 取模)
        device_id = (worker_id - 1) % 4
        device = torch.device(f"cuda:{device_id}")
        
        file_name = os.path.basename(file_path)
        out_file_name = file_name.replace('forward_inputs_', 'open_lands_outputs_')
        if out_file_name == file_name:
            out_file_name = 'open_lands_' + file_name
        out_file_path = os.path.join(output_dir, out_file_name)

        # ==========================================
        # 1. 检查文件是否已存在 (断点续传)
        # ==========================================
        if os.path.exists(out_file_path):
            print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {out_file_name} already exists. Skipping...")
            return True, file_path, f"Skipped: {out_file_name} exists"

        print(f"\n[{time.strftime('%H:%M:%S')}] [File {file_idx}/{total_files}] Loading {file_name} on GPU {device_id}...")
        t_file_start = time.time()
        
        # 初始化模型 (多进程环境下必须在子进程中独立初始化)
        rtm = SoilRTM_TB_Calculator(def_da_rtm_rough=2).to(device)
        rtm.eval()
        dz_soi = rtm.dz_soi.to(device)
        
        wtot = 0.0175 + 0.0276  # 表面前两层总厚度 (0.0451 m)
        denh2o = 1000.0         # 水密度
        denice = 917.0          # 冰密度

        # 恢复老代码的大吞吐分块策略
        time_chunk_size = 8784  
        batch_size = 15000

        # ==========================================
        # 2. 读取原始数据并进行内部筛选
        # ==========================================
        ds_in = nc.Dataset(file_path, 'r')
        ds_in.set_auto_mask(False)
        
        open_lands_classes = {6, 7, 8, 9, 10, 11, 12, 14, 16}
        patchclass_all = ds_in.variables['patchclass'][:]
        mask = np.isin(patchclass_all, list(open_lands_classes))
        
        if not np.any(mask):
            ds_in.close()
            return True, file_path, f"No valid open lands patches in {file_name}"

        local_p_indices = np.where(mask)[0]
        n_patches = len(local_p_indices)
        num_times = ds_in.variables['time'].shape[0]
        
        # 🚀 修复：先完整读入内存，再利用 local_p_indices 进行切片
        wf_clay_raw = ds_in.variables['wf_clay'][:]
        wf_clay_raw = wf_clay_raw[local_p_indices, :]
        
        wf_sand_raw = ds_in.variables['wf_sand'][:]
        wf_sand_raw = wf_sand_raw[local_p_indices, :]
        
        porsl_raw   = ds_in.variables['porsl'][:]
        porsl_raw   = porsl_raw[local_p_indices, :]
        
        pclass_raw  = patchclass_all[local_p_indices]
        
        lon_raw     = ds_in.variables['lon'][:]
        lon_raw     = lon_raw[local_p_indices]
        
        lat_raw     = ds_in.variables['lat'][:]
        lat_raw     = lat_raw[local_p_indices]
        
        if 'sat_theta' in ds_in.variables:
            theta_raw = ds_in.variables['sat_theta'][:]
            theta_raw = theta_raw[local_p_indices]
        else:
            theta_raw = np.full(n_patches, 40.0 * math.pi / 180.0)
            
        if 'sat_fghz' in ds_in.variables:
            fghz_raw = ds_in.variables['sat_fghz'][:]
            fghz_raw = fghz_raw[local_p_indices]
        else:
            fghz_raw = np.full(n_patches, 1.4)
        
        # 创建独立的输出 NC 文件
        ds_out = nc.Dataset(out_file_path, 'w', format='NETCDF4')
        ds_out.createDimension('time', num_times)
        ds_out.createDimension('patch', n_patches)
        
        v_time = ds_out.createVariable('time', 'f8', ('time',))
        ds_out.createVariable('patch_idx', 'i4', ('patch',))[:] = local_p_indices
        ds_out.createVariable('lon', 'f8', ('patch',))[:] = lon_raw
        ds_out.createVariable('lat', 'f8', ('patch',))[:] = lat_raw
        ds_out.createVariable('patchclass', 'i4', ('patch',))[:] = pclass_raw
        
        v_time[:] = ds_in.variables['time'][:]

        schemes = ['wilheit', 'holmes', 'wigneron', 'shuyueliu', 'landem']
        pols = ['h', 'v']
        out_vars = {}
        for s in schemes:
            for p in pols:
                v_name = f'tb_{s}_{p}'
                out_vars[v_name] = ds_out.createVariable(v_name, 'f4', ('time', 'patch'), zlib=True)

        # 计算表面层与深层静态属性的加权
        wf_clay_surf = (wf_clay_raw[:, 0] * 0.0175 + wf_clay_raw[:, 1] * 0.0276) / wtot
        wf_sand_surf = (wf_sand_raw[:, 0] * 0.0175 + wf_sand_raw[:, 1] * 0.0276) / wtot
        porsl_surf   = (porsl_raw[:, 0] * 0.0175 + porsl_raw[:, 1] * 0.0276) / wtot

        base_clay_surf = torch.tensor(wf_clay_surf, dtype=torch.float32, device=device)
        base_sand_surf = torch.tensor(wf_sand_surf, dtype=torch.float32, device=device)
        base_porsl_surf = torch.tensor(porsl_surf, dtype=torch.float32, device=device)
        base_pclass = torch.tensor(pclass_raw, dtype=torch.long, device=device)
        base_theta = torch.tensor(theta_raw, dtype=torch.float32, device=device)
        base_fghz = torch.tensor(fghz_raw, dtype=torch.float32, device=device)

        # ==========================================
        # 3. GPU 推理核心循环 (外层分块读取，内层显存切分)
        # ==========================================
        with torch.no_grad():
            for t_start in range(0, num_times, time_chunk_size):
                t_end = min(t_start + time_chunk_size, num_times)
                cur_T = t_end - t_start
                
                # 🚀 修复：先全量读取时间块的大矩阵，避免底层随机碎片寻道
                t_soisno_raw = ds_in.variables['t_soisno'][t_start:t_end, :, 5:15]
                wliq_soisno_raw = ds_in.variables['wliq_soisno'][t_start:t_end, :, 5:15]
                wice_soisno_raw = ds_in.variables['wice_soisno'][t_start:t_end, :, 5:15]

                # 在内存中进行高速掩码切片
                t_soisno_chunk = t_soisno_raw[:, local_p_indices, :]
                wliq_soisno_chunk = wliq_soisno_raw[:, local_p_indices, :]
                wice_soisno_chunk = wice_soisno_raw[:, local_p_indices, :]

                # 仅在 CPU 端创建 Full Tensor，不在此处放入 device
                t_soi_full = torch.tensor(t_soisno_chunk, dtype=torch.float32).reshape(-1, 10)
                wliq_soi_full = torch.tensor(wliq_soisno_chunk, dtype=torch.float32).reshape(-1, 10)
                wice_soi_full = torch.tensor(wice_soisno_chunk, dtype=torch.float32).reshape(-1, 10)
                
                # 及时释放内存
                del t_soisno_raw, wliq_soisno_raw, wice_soisno_raw
                del t_soisno_chunk, wliq_soisno_chunk, wice_soisno_chunk
                
                total_samples = cur_T * n_patches
                chunk_res = {s: [] for s in schemes}

                for i in range(0, total_samples, batch_size):
                    end_idx = min(i + batch_size, total_samples)
                    
                    # 仅将当前 Batch 的数据推入对应的 GPU，non_blocking 加速传输
                    t_soi = t_soi_full[i:end_idx].to(device, non_blocking=True)
                    wliq = wliq_soi_full[i:end_idx].to(device, non_blocking=True)
                    wice = wice_soi_full[i:end_idx].to(device, non_blocking=True)
                    
                    patch_match_indices = torch.arange(i, end_idx, device=device) % n_patches
                    
                    b_clay = base_clay_surf[patch_match_indices]
                    b_sand = base_sand_surf[patch_match_indices]
                    b_porsl = base_porsl_surf[patch_match_indices]
                    b_pclass = base_pclass[patch_match_indices]
                    b_theta = base_theta[patch_match_indices]
                    b_fghz = base_fghz[patch_match_indices]
                    
                    # ------------------- 物理量加权计算 (对齐 Fortran) -------------------
                    liq_soi_profile = wliq / (dz_soi * denh2o)
                    ice_soi_profile = wice / (dz_soi * denice)
                    
                    total_water_profile = liq_soi_profile + ice_soi_profile
                    ffrz_profile = torch.where(total_water_profile > 0, 
                                               ice_soi_profile / total_water_profile, 
                                               torch.zeros_like(total_water_profile))
                    
                    liq_surf = (wliq[:, 0] + wliq[:, 1]) / (wtot * denh2o)
                    ice_surf = (wice[:, 0] + wice[:, 1]) / (wtot * denice)
                    
                    total_surf = liq_surf + ice_surf
                    ffrz_surf = torch.where(total_surf > 0, 
                                            ice_surf / total_surf, 
                                            torch.zeros_like(total_surf))
                    
                    t_surf_k = (t_soi[:, 0] * 0.0175 + t_soi[:, 1] * 0.0276) / wtot
                    t_deep_k = (t_soi[:, 6] * (0.8289 - 0.5) + t_soi[:, 7] * (1.0 - 0.8289)) / 0.5
                    
                    b_is_desert = (liq_surf < 0.02) & (b_sand > 90.0)
                    # ---------------------------------------------------------------------

                    # 执行张量辐射传输算子 (加入细粒度日志打印与防崩溃捕获)
                    try:
                        tb_w = rtm.calculate_tb_soil_CMEM('wilheit', t_surf_k, t_deep_k, t_soi, liq_surf, liq_soi_profile, b_clay, ffrz_surf, ffrz_profile, b_theta, b_fghz, b_is_desert, b_pclass)
                    except Exception as e:
                        print(f"❌ [Batch {i}:{end_idx}] Error in Wilheit scheme: {e}")
                        raise

                    try:
                        tb_h = rtm.calculate_tb_soil_CMEM('holmes', t_surf_k, t_deep_k, t_soi, liq_surf, liq_soi_profile, b_clay, ffrz_surf, ffrz_profile, b_theta, b_fghz, b_is_desert, b_pclass)
                    except Exception as e:
                        print(f"❌ [Batch {i}:{end_idx}] Error in Holmes scheme: {e}")
                        raise

                    try:
                        tb_wg = rtm.calculate_tb_soil_CMEM('wigneron', t_surf_k, t_deep_k, t_soi, liq_surf, liq_soi_profile, b_clay, ffrz_surf, ffrz_profile, b_theta, b_fghz, b_is_desert, b_pclass)
                    except Exception as e:
                        print(f"❌ [Batch {i}:{end_idx}] Error in Wigneron scheme: {e}")
                        raise

                    try:
                        tb_s = rtm.calculate_tb_soil_shuyueliu(t_soi, liq_soi_profile, b_clay, ffrz_profile, b_theta, b_fghz, b_is_desert, b_pclass, b_sand, b_porsl)
                    except Exception as e:
                        print(f"❌ [Batch {i}:{end_idx}] Error in ShuyueLiu scheme: {e}")
                        raise

                    try:
                        tb_l = rtm.calculate_tb_soil_LandEM(t_surf_k, t_deep_k, liq_surf, b_clay, ffrz_surf, b_theta, b_fghz, b_is_desert, b_pclass)
                    except Exception as e:
                        print(f"❌ [Batch {i}:{end_idx}] Error in LandEM scheme: {e}")
                        raise
                    
                    chunk_res['wilheit'].append(tb_w.cpu())
                    chunk_res['holmes'].append(tb_h.cpu())
                    chunk_res['wigneron'].append(tb_wg.cpu())
                    chunk_res['shuyueliu'].append(tb_s.cpu())
                    chunk_res['landem'].append(tb_l.cpu())

                # 将当前时间块的数据拼接并写入 NetCDF
                for s in schemes:
                    concat_tb = torch.cat(chunk_res[s], dim=1).numpy()
                    reshaped_tb = concat_tb.reshape(2, cur_T, n_patches)
                    out_vars[f'tb_{s}_h'][t_start:t_end, :] = reshaped_tb[0]
                    out_vars[f'tb_{s}_v'][t_start:t_end, :] = reshaped_tb[1]

        ds_in.close()
        ds_out.close()
        print(f"✅ [GPU {device_id} | {time.strftime('%H:%M:%S')}] Finished {out_file_name} | Cost: {time.time() - t_file_start:.2f} s")
        return True, file_path, f"Processed {file_name} -> GPU {device_id}"

    except Exception as e:
        print(f"❌ [GPU {device_id}] Error processing {file_path}: {e}")
        return False, file_path, f"Error in {os.path.basename(file_path)}: {e}"
        
    finally:
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_file.close()

if __name__ == '__main__':
    # 强制在子进程中重新分配 CUDA 内存池以防止死锁
    mp.set_start_method('spawn', force=True)
    
    # --- 1. 路径与配置 ---
    nc_dir = '/home/liusy/CoLM/outputs/global_veg_wigneron/forward_inputs_folder'
    output_dir = '/home/liusy/research_lists/2026-07-15_research_list/compare_cmem_landem_ours/open_lands_optimized_results'
    os.makedirs(output_dir, exist_ok=True)
    
    # 获取所有的 NC 文件列表
    nc_files = sorted(glob.glob(os.path.join(nc_dir, 'forward_inputs_worker*.nc')))
    total_files = len(nc_files)
    
    MAX_WORKERS = 4  # 4进程对应 4张 5090
    
    print(f"\n=======================================================")
    print(f"🚀 开始多卡高并发提取，共发现 {total_files} 个 NetCDF 文件")
    print(f"分配策略: {MAX_WORKERS} 个并发 Worker 均匀分布在 4 张 5090 显卡上")
    print(f"日志状态: 进程内计算细节将被静默记录到 `log_files_compare` 目录")
    print(f"=======================================================\n")

    t_global_start = time.time()

    # 启用 ProcessPoolExecutor 强力调度
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # 将任务提交入池
        futures = {
            executor.submit(worker_process, idx, file_path, total_files, output_dir): file_path 
            for idx, file_path in enumerate(nc_files, 1)
        }
        
        # tqdm 在主线程接管进度，展示极其整洁的 UI
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="🚀 全局计算进度"):
            success, file_path, msg = future.result()
            
            # 只有遇到错误时，才打破进度条在主终端报错提醒你
            if not success:
                tqdm.write(f"❌ 警告: {msg}")

    print(f"\n🎉 All files processed successfully! Total Elapsed Time: {(time.time() - t_global_start) / 60:.2f} mins.")