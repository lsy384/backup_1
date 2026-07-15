import torch
import torch.nn as nn
import math
import numpy as np
import pands as pd


def eff_soil_temp_Wilheit(dz_soi, t_soi, eps, theta, lamcm):
        """
        修正后的 Wilheit 1978 相干模型。
        在张量化操作中严格复现了 Fortran 物理衰减截断 (NMAX) 的逻辑，避免数值爆炸。
        """
        batch_size = eps.shape[0]
        nlay = eps.shape[1]
        N = nlay + 1
        device = eps.device

        # --- 数据精度建议 ---
        # 注: Fortran原版使用了 JPRB (通常是float64/complex128)。如果复现精度要求高，
        # 建议将 dtype 提升为 torch.complex128 和 torch.float64。
        dtype_c = torch.complex64 
        dtype_f = torch.float32

        # 1. 计算折射率 (大气层折射率为 1.0)
        CN = torch.zeros((batch_size, N), dtype=dtype_c, device=device)
        CN[:, 0] = 1.0 + 0.0j
        CN[:, 1:] = torch.sqrt(eps)

        DEL = torch.zeros((batch_size, N), dtype=dtype_f, device=device)
        DEL[:, 1:] = dz_soi * 100.0  # 假设 dz_soi 是米，转为 cm [cite: 27]

        # 2. 计算界面传播矩阵 (CP) 及其动态截断点 NMAX
        S = torch.sin(theta)
        CP = torch.ones((batch_size, N), dtype=dtype_c, device=device)
        
        # NMAX_idx 用于追踪每个 batch 样本的电磁波有效穿透层
        NMAX_idx = torch.full((batch_size,), N-1, dtype=torch.long, device=device)
        active_mask = torch.ones(batch_size, dtype=torch.bool, device=device)

        for i in range(1, N-1):
            CS = CN[:, 0] * S / CN[:, i]
            CC = torch.sqrt(1.0 + 0.0j - CS*CS)
            ARG = DEL[:, i] * 2.0 * math.pi / lamcm
            CARG = 2.0 * ARG * CN[:, i] * CC * 1.0j
            
            CP_next = torch.exp(CARG) * CP[:, i-1]
            CP[:, i] = torch.where(active_mask, CP_next, CP[:, i-1])
            
            # 严格对齐 Fortran: IF (ABS(CP(I)) < 0.0001) EXIT
            below_thresh = torch.abs(CP[:, i]) < 0.0001
            crossed_now = active_mask & below_thresh
            
            # 记录刚刚穿过阈值的层索引
            NMAX_idx = torch.where(crossed_now, torch.tensor(i + 1, device=device), NMAX_idx)
            active_mask = active_mask & (~below_thresh)

        # 3. 逆向计算电场 (Backward loop)
        CEP_h = torch.zeros((batch_size, N), dtype=dtype_c, device=device)
        CEM_h = torch.zeros((batch_size, N), dtype=dtype_c, device=device)
        CEP_v = torch.zeros((batch_size, N), dtype=dtype_c, device=device)
        CEM_v = torch.zeros((batch_size, N), dtype=dtype_c, device=device)

        # Fortran: 仅在 NMAX 处初始化边界条件为 (1, 0)
        batch_indices = torch.arange(batch_size, device=device)
        CEP_h[batch_indices, NMAX_idx] = 1.0 + 0.0j
        CEP_v[batch_indices, NMAX_idx] = 1.0 + 0.0j

        # 逆向传播，从最大可能层往回算
        for j in range(N-2, -1, -1):
            active_j = j < NMAX_idx  # 只计算处于有效穿透深度的物理层
            
            # 为无效层提供一个安全的除数兜底，防止 PyTorch 后台报错，但这部分数据会被屏蔽
            safe_cp = torch.where(active_j, CP[:, j], torch.tensor(1.0, dtype=dtype_c, device=device))

            CSJ = CN[:, 0] * S / CN[:, j]
            CCJ = torch.sqrt(1.0 + 0.0j - CSJ*CSJ)
            CSJP1 = CN[:, 0] * S / CN[:, j+1]
            CCJP1 = torch.sqrt(1.0 + 0.0j - CSJP1*CSJP1)

            # --- 水平极化 (HPLD) ---
            CA_h = 2.0 * CN[:, j] * CCJ / (CN[:, j]*CCJ + CN[:, j+1]*CCJP1)
            CB_h = (CN[:, j]*CCJ - CN[:, j+1]*CCJP1) / ((CN[:, j]*CCJ + CN[:, j+1]*CCJP1) * safe_cp)
            
            CEP_h_j = CEP_h[:, j+1]/CA_h + CB_h*CEM_h[:, j+1]/CA_h
            CEM_h_j = CEM_h[:, j+1] + (CEP_h[:, j+1] - CEP_h_j)*safe_cp

            # 利用 active_j 门控：如果是无效层，保持电场为 0
            CEP_h[:, j] = torch.where(active_j, CEP_h_j, CEP_h[:, j])
            CEM_h[:, j] = torch.where(active_j, CEM_h_j, CEM_h[:, j])

            # --- 垂直极化 (VPLD) ---
            CD_v = 2.0 * CN[:, j] * CCJ
            CA_v = CN[:, j]*CCJP1 + CN[:, j+1]*CCJ
            CB_v = CN[:, j]*CCJP1 - CN[:, j+1]*CCJ
            
            CEP_v_j = CA_v*CEP_v[:, j+1]/CD_v + CB_v*CEM_v[:, j+1]/(CD_v*safe_cp)
            CR_v = CN[:, j+1] / CN[:, j]
            CEM_v_j = CR_v*CEM_v[:, j+1] + (CEP_v_j - CEP_v[:, j+1]*CR_v)*safe_cp
            
            CEP_v[:, j] = torch.where(active_j, CEP_v_j, CEP_v[:, j])
            CEM_v[:, j] = torch.where(active_j, CEM_v_j, CEM_v[:, j])

        # 归一化电场
        CX_h = CEP_h[:, 0]
        CX_v = CEP_v[:, 0]
        safe_CX_h = torch.where(torch.abs(CX_h) > 0, CX_h, torch.tensor(1.0, dtype=dtype_c, device=device))
        safe_CX_v = torch.where(torch.abs(CX_v) > 0, CX_v, torch.tensor(1.0, dtype=dtype_c, device=device))
        
        for j in range(N):
            CEP_h[:, j] = CEP_h[:, j] / safe_CX_h
            CEM_h[:, j] = CEM_h[:, j] / safe_CX_h
            CEP_v[:, j] = CEP_v[:, j] / safe_CX_v
            CEM_v[:, j] = CEM_v[:, j] / safe_CX_v

        # 4. 计算各层吸收率 (fa)
        fa_h = torch.zeros((batch_size, nlay), dtype=dtype_f, device=device)
        fa_v = torch.zeros((batch_size, nlay), dtype=dtype_f, device=device)
        cos_theta = torch.cos(theta)

        # 严格对齐 Fortran: 对于 NMAX 及以下的层，将 CP 重置为极小虚数 
        j_indices = torch.arange(N, device=device).unsqueeze(0)
        mask_nmax = j_indices >= NMAX_idx.unsqueeze(1)
        CP = torch.where(mask_nmax, torch.tensor(1e-15j, dtype=dtype_c, device=device), CP)

        for j in range(1, N):
            active_j = j <= NMAX_idx  # 吸收率依然受 NMAX 截断控制
            
            CS = S / CN[:, j]
            CC = torch.sqrt(1.0 + 0.0j - CS*CS)
            
            R = torch.abs(CP[:, j])
            S_abs = torch.abs(CP[:, j-1])
            
            # 数值保护，因为有倒数项
            safe_R = torch.where(R > 1e-25, R, torch.tensor(1e-25, device=device))
            safe_S = torch.where(S_abs > 1e-25, S_abs, torch.tensor(1e-25, device=device))

            # --- HPLD 吸收 ---
            E2_h = (safe_S - safe_R) * torch.abs(CEP_h[:, j])**2 + (1.0/safe_R - 1.0/safe_S) * torch.abs(CEM_h[:, j])**2
            DP_h = E2_h * (CN[:, j]*CC).real / cos_theta
            CXP_h = CEP_h[:, j] * torch.conj(CEM_h[:, j])
            X_h = 2.0 * ((CN[:, j]*CC / cos_theta).imag) * ( (CXP_h * CP[:, j-1]/safe_S).imag - (CXP_h * CP[:, j]/safe_R).imag )
            
            fa_h_val = DP_h - X_h
            fa_h[:, j-1] = torch.where(active_j, fa_h_val.real, torch.zeros_like(fa_h_val.real))
            
            # --- VPLD 吸收 ---
            E2_v = (safe_S - safe_R) * torch.abs(CEP_v[:, j])**2 + (1.0/safe_R - 1.0/safe_S) * torch.abs(CEM_v[:, j])**2
            DP_v = E2_v * (CN[:, j]*CC).real / cos_theta
            CXP_v = CEP_v[:, j] * torch.conj(CEM_v[:, j])
            X_v = 2.0 * ((CN[:, j]*CC / cos_theta).imag) * ( (CXP_v * CP[:, j-1]/safe_S).imag - (CXP_v * CP[:, j]/safe_R).imag )
            
            fa_v_val = DP_v - X_v
            fa_v[:, j-1] = torch.where(active_j, fa_v_val.real, torch.zeros_like(fa_v_val.real))

        # 5. 计算有效温度和表面反射率
        sum_fa_h = torch.sum(fa_h, dim=1)
        teff_h = torch.where(sum_fa_h == 0, t_soi[:, 0], torch.sum(fa_h * t_soi, dim=1) / sum_fa_h)
        
        sum_fa_v = torch.sum(fa_v, dim=1)
        teff_v = torch.where(sum_fa_v == 0, t_soi[:, 0], torch.sum(fa_v * t_soi, dim=1) / sum_fa_v)

        # r_s(1) = real(cp(1)) 中的 R = ABS(CEM(1))**2. * REAL(CN(1)) [cite: 46]
        r_h = torch.abs(CEM_h[:, 0])**2 * CN[:, 0].real
        r_v = torch.abs(CEM_v[:, 0])**2 * CN[:, 0].real
        
        return r_h, r_v, teff_h, teff_v    
    
    
def eff_soil_temp_Wigneron2001(wc_surf, t_surf, t_deep):
        """
        Wigneron 2001 scheme
        """
        w0 = 0.30
        bw = 0.30
        
        # 🌟 防爆修复 3：在除以 w0 和求 bw 次幂之前，严格限制 wc_surf 为正数
        wc_safe = torch.clamp(wc_surf, min=1e-6)
        C = torch.clamp((wc_safe / w0)**bw, min=0.001)
        teff = t_deep + (t_surf - t_deep) * C
        return teff
    
    
def eff_soil_temp_Holmes2006(eps_surf, t_surf, t_deep, return_C=False):
        """
        Holmes 2006 model (基于介电常数比值的参数化)
        eps_surf: 表层复介电常数 (complex)
        """
        eps_r = eps_surf.real
        eps_i = torch.abs(eps_surf.imag)
        
        # 2003-2004 interannual calibration parameters
        b_param = 0.05      # 0.87
        eps0_param = 0.08

        # C(eps) = ((eps'' / eps') / eps0_param)^b
        eps_ratio = eps_i / eps_r
        C = (eps_ratio / eps0_param) ** b_param
        C = torch.clamp(C, min=0.001)
        
        teff = t_deep + (t_surf - t_deep) * C
        if return_C:
            return teff, C
        return teff        
    
