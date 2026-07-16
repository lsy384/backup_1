import torch
import torch.nn as nn
import math

class SoilRTM_TB_Calculator(nn.Module):
    def __init__(self, def_da_rtm_rough=2):
        """
        初始化土壤亮温计算器
        :param def_da_rtm_rough: 粗糙度模型选项 (默认为 2, 对应使用 hr_SMAP 且 nrh=2.0, nrv=2.0)
        """
        super(SoilRTM_TB_Calculator, self).__init__()
        
        # --- 物理常数与基本配置 ---
        self.pi = math.pi
        self.tfrz = 273.15           # freezing temperature [K]
        self.C = 299792458.0         # speed of light in vacuum [m/s]
        self.eps_0 = 8.854187817e-12 # vacuum permittivity [F/m]
        self.eps_w_inf = 4.9         # high-frequency limit of water permittivity

        # --- 默认土壤各层厚度 (m) ---
        self.register_buffer('dz_soi', torch.tensor([0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 
                                                     0.2038, 0.3360, 0.5539, 0.9133, 1.5058]))
        
        # --- 选项配置与常量 ---
        self.def_da_rtm_rough = def_da_rtm_rough  
        self.rgh_surf = 2.2          # 默认表面粗糙度参数
        
        # --- 注册粗糙度经验参数表 (依附于模型，支持设备迁移) ---
        self.register_buffer('hr_SMAP', torch.tensor([0.0, 
            0.160, 0.160, 0.160, 0.160, 0.160, 0.110, 0.110, 0.125, 0.156, 0.156, 
            0.100, 0.108, 0.000, 0.130, 0.000, 0.150, 0.000]))
        self.register_buffer('hr_SMOS', torch.tensor([0.0, 
            0.300, 0.300, 0.300, 0.300, 0.300, 0.100, 0.100, 0.100, 0.100, 0.100, 
            0.100, 0.100, 0.100, 0.100, 0.000, 0.100, 0.000]))
        self.register_buffer('hr_P16', torch.tensor([0.0, 
            0.350, 0.460, 0.430, 0.450, 0.410, 0.260, 0.170, 0.350, 0.230, 0.130, 
            0.020, 0.170, 0.190, 0.220, 0.000, 0.020, 0.000]))
        
        # --- 初始化中间状态变量 ---
        self.t_eff = None
        self.t_eff_wilheit = None
        self.t_eff_holmes = None
        self.t_eff_wigneron = None
        self.r_s = None
        self.r_r = None 
        self.tb_soil = None
        self.tb_soil_d = None
        self.tb_soil_nd = None

    def diel_soil_M09(self, wc, t, wf_clay, f):
        """
        Mironov 2009 土壤介电常数模型 (仅使用粘土含量)
        wc      : 体积含水量 (液态) [m3/m3]，任意形状张量
        t       : 土壤温度 [K]，与 wc 形状相同
        wf_clay : 粘土质量百分含量 [%]，需要可广播到 wc 形状
        f       : 频率 [Hz]，需要可广播到 wc 形状
        返回    : 复介电常数张量，与 wc 形状相同
        """
        device = wc.device
        dtype = wc.dtype

        wf_clay_frac = wf_clay / 100.0
        tk = t + self.tfrz

        nd = 1.634 - 0.539 * wf_clay_frac + 0.2748 * (wf_clay_frac ** 2)        
        kd = 0.03952 - 0.04038 * wf_clay_frac                                   

        mvt = 0.02863 + 0.30673 * wf_clay_frac                                  
        ts = 20.0
        
        e0b = 79.8 - 85.4 * wf_clay + 32.7 * (wf_clay_frac ** 2)                
        Bb = 8.67e-19 - 0.00126 * wf_clay_frac + 0.00184 * (wf_clay_frac ** 2) \
             - 9.77e-10 * (wf_clay ** 3) - 1.39e-15 * (wf_clay ** 4)            
        Bsgb = 0.0028 + 0.02094e-2 * wf_clay - 0.01229e-4 * (wf_clay ** 2) \
               - 5.03e-22 * (wf_clay ** 3) + 4.163e-24 * (wf_clay ** 4)         

        # 确保参与 torch.log 的基础项是 Tensor
        e0b_tensor = torch.as_tensor(e0b, device=device, dtype=dtype)
        Fb = torch.log((e0b_tensor - 1.0) / (e0b_tensor + 2.0))                               
        exp_term_b = torch.exp(Fb - Bb * (t - ts))
        eb0 = (1.0 + 2.0 * exp_term_b) / (1.0 - exp_term_b)                     

        dHbR = 1467.0 + 2697e-2 * wf_clay - 980e-4 * (wf_clay ** 2) \
               + 1.368e-10 * (wf_clay ** 3) - 8.61e-13 * (wf_clay ** 4)         
        dSbR = 0.888 + 9.7e-2 * wf_clay - 4.262e-4 * (wf_clay ** 2) \
               + 6.79e-21 * (wf_clay ** 3) + 4.263e-22 * (wf_clay ** 4)         
               
        taub = 48e-12 * torch.exp(dHbR / tk - dSbR) / tk                        

        sigmabt = 0.3112 + 0.467e-2 * wf_clay                                   
        sigmab = sigmabt + Bsgb * (t - ts)                                      

        # 将 e0u 挂在 Tensor 上，防止后续 torch.log 报 TypeError
        e0u = torch.tensor(100.0, device=device, dtype=dtype)                                                             
        Bu = 1.11e-4 - 1.603e-7 * wf_clay + 1.239e-9 * (wf_clay ** 2) \
             + 8.33e-13 * (wf_clay ** 3) - 1.007e-14 * (wf_clay ** 4)           
        Bsgu = 0.00108 + 0.1413e-2 * wf_clay - 0.2555e-4 * (wf_clay ** 2) \
               + 0.2147e-6 * (wf_clay ** 3) - 0.0711e-8 * (wf_clay ** 4)        

        Fu = torch.log((e0u - 1.0) / (e0u + 2.0))                               
        exp_term_u = torch.exp(Fu - Bu * (t - ts))
        eu0 = (1.0 + 2.0 * exp_term_u) / (1.0 - exp_term_u)                     

        dHuR = 2231.0 - 143.1e-2 * wf_clay + 223.2e-4 * (wf_clay ** 2) \
               - 142.1e-6 * (wf_clay ** 3) + 27.14e-8 * (wf_clay ** 4)          
        dSuR = 3.649 - 0.4894e-2 * wf_clay + 0.763e-4 * (wf_clay ** 2) \
               - 0.4859e-6 * (wf_clay ** 3) + 0.0928e-8 * (wf_clay ** 4)        
               
        tauu = 48e-12 * torch.exp(dHuR / tk - dSuR) / tk                        

        sigmaut = 0.05 + 1.4 * (1.0 - (1.0 - wf_clay * 1e-2) ** 4.664)          
        sigmau = sigmaut + Bsgu * (t - ts)                                      

        cxb = (eb0 - self.eps_w_inf) / (1.0 + (2.0 * self.pi * f * taub)**2)    
        eb_r = self.eps_w_inf + cxb                                             
        eb_i = cxb * (2.0 * self.pi * f * taub) + sigmab / (2.0 * self.pi * self.eps_0 * f) 

        cxu = (eu0 - self.eps_w_inf) / (1.0 + (2.0 * self.pi * f * tauu)**2)    
        eu_r = self.eps_w_inf + cxu                                             
        eu_i = cxu * (2.0 * self.pi * f * tauu) + sigmau / (2.0 * self.pi * self.eps_0 * f) 

        # 挂在 Tensor 上的开方操作
        sqrt_2 = torch.sqrt(torch.tensor(2.0, device=device, dtype=dtype))
        nb = torch.sqrt(torch.sqrt(eb_r**2 + eb_i**2) + eb_r) / sqrt_2          
        kb = torch.sqrt(torch.sqrt(eb_r**2 + eb_i**2) - eb_r) / sqrt_2          
        nu = torch.sqrt(torch.sqrt(eu_r**2 + eu_i**2) + eu_r) / sqrt_2          
        ku = torch.sqrt(torch.sqrt(eu_r**2 + eu_i**2) - eu_r) / sqrt_2          

        is_le_mvt = wc <= mvt
        nm_le = nd + (nb - 1.0) * wc                                             
        km_le = kd + kb * wc                                                     
        
        nm_gt = nd + (nb - 1.0) * mvt + (nu - 1.0) * (wc - mvt)                 
        km_gt = kd + kb * mvt + ku * (wc - mvt)                                 
        
        nm = torch.where(is_le_mvt, nm_le, nm_gt)
        km = torch.where(is_le_mvt, km_le, km_gt)

        eps_r = nm**2 - km**2                                                    
        eps_i = 2.0 * nm * km                                                    

        return torch.complex(eps_r, eps_i)

    def calculate_tb_soil_CMEM(self, scheme, t_surf_k, t_deep_k, t_soi_k,
                          liq_surf, liq_soi_profile, wf_clay_surf,
                          ffrz_surf, ffrz_profile,
                          theta, fghz, is_desert, patchclass):
        """
        主函数：根据选择的有效温度方案计算土壤亮温 tb_soil
        
        新增参数:
        liq_surf         : 表面液态水体积含水量 [m3/m3]  (batch,)
        liq_soi_profile  : 各层液态水体积含水量 [m3/m3]  (batch, 10)
        wf_clay_surf     : 表面粘土质量百分含量 [%]      (batch,)
        ffrz_surf        : 表面冻结比例 [0~1]           (batch,)
        ffrz_profile     : 各层冻结比例 [0~1]           (batch, 10)
        """
        # 1. 计算物理参数
        f = fghz * 1e9                 
        lam = self.C / f               
        k = 2 * self.pi / lam          
        kcm = k / 100.0                
        kr = k * (0.5 * 1e-3)          
        lamcm = lam * 100.0

        # 冻结冰的介电常数 (M09 仅计算液态，冻土混合在此处理)
        eps_f = torch.complex(torch.full_like(liq_surf, 5.0), torch.full_like(liq_surf, 0.5))

        # 2. 统一使用 M09 计算表面及剖面介电常数
        eps_soil_surf_liq = self.diel_soil_M09(liq_surf, t_surf_k, wf_clay_surf, f)
        eps_soil_surf = eps_soil_surf_liq * (1.0 - ffrz_surf) + eps_f * ffrz_surf
        
        # 3. 根据方案计算有效温度 T_eff，并保存至对应的 self 变量
        if scheme.lower() == 'wilheit':
            r_h, r_v, teff_h, teff_v = self.eff_soil_temp_Wilheit(
                t_soi_k, liq_soi_profile, wf_clay_surf, f, ffrz_profile, theta, lamcm)
            self.t_eff_wilheit = torch.stack([teff_h, teff_v])
            self.t_eff = self.t_eff_wilheit
            
        elif scheme.lower() == 'holmes':
            teff = self.eff_soil_temp_Holmes2006(eps_soil_surf, t_surf_k, t_deep_k)
            self.t_eff_holmes = torch.stack([teff, teff])
            self.t_eff = self.t_eff_holmes
            
        elif scheme.lower() == 'wigneron':
            teff = self.eff_soil_temp_Wigneron2001(liq_surf, t_surf_k, t_deep_k)
            self.t_eff_wigneron = torch.stack([teff, teff])
            self.t_eff = self.t_eff_wigneron
            
        else:
            raise ValueError("不支持的有效温度方案。请选择 'wilheit', 'holmes', 或 'wigneron'。")

        # 4. 计算菲涅尔平滑表面反射率 (Fresnel reflectivity)
        g = torch.sqrt(eps_soil_surf - torch.sin(theta)**2)
        r_s_h = torch.abs((torch.cos(theta) - g)/(torch.cos(theta) + g))**2
        r_s_v = torch.abs((torch.cos(theta)*eps_soil_surf - g)/(torch.cos(theta)*eps_soil_surf + g))**2
        
        self.r_s = torch.stack([r_s_h, r_s_v])

        # 5. 计算粗糙表面反射率 (Rough reflectivity)，移除了 hr_pred
        self.r_r = self.rough_reflectivity(is_desert, patchclass, self.r_s, theta, fghz, kcm, self.rgh_surf) 

        # 6. 计算最终土壤亮温 (分荒漠和非荒漠)
        self.tb_soil_d = self.desert(self.t_eff, self.r_r, eps_soil_surf, kr)
        self.tb_soil_nd = self.t_eff * (1 - self.r_r)
        
        is_desert_2d = is_desert.unsqueeze(0).expand(2, -1)
        self.tb_soil = torch.where(is_desert_2d, self.tb_soil_d, self.tb_soil_nd)
        
        return self.tb_soil
  
    def calculate_tb_soil_shuyueliu(self, t_soi_k, liq_soi_profile, wf_clay_surf, ffrz_profile,
                                    theta, fghz, is_desert, patchclass, wf_sand, porsl):
        """
        Shuyue_LIU 离散多层 RTE 土壤亮温计算方法。
        
        该方案考虑了土壤剖面体散射的双流近似辐射传输模型，并使用矩阵法严格求解。
        输入接口已与 CMEM 方案对齐，直接输入土壤廓线和地表基础变量，二级物理量（介电常数及粗糙地表反射率）在内部动态计算。
        依据要求，边界条件排除了大气下行辐射，输出仅包含“发射率 * 有效温度”的土壤自身向上发射项。

        参数说明:
        t_soi_k          : [tensor] 土壤各层物理温度，形状 (batch, n_layers)，单位: K
        liq_soi_profile  : [tensor] 各层液态水体积含水量，形状 (batch, n_layers)，单位: m3/m3
        wf_clay_surf     : [tensor] 表面粘土质量百分含量，形状 (batch,)，单位: %
        ffrz_profile     : [tensor] 各层冻结比例，形状 (batch, n_layers)，范围: 0~1
        theta            : [tensor/float] 观测天顶角，单位: rad
        fghz             : [tensor/float] 微波频率，单位: GHz
        is_desert        : [tensor] 是否为荒漠的布尔掩码，形状 (batch,)
        patchclass       : [tensor] 地表覆盖类型索引，形状 (batch,)
        wf_sand          : [tensor/float] 沙粒质量分数
        porsl            : [tensor/float] 土壤孔隙度

        返回:
        tb_soil          : [tensor] 形状为 (2, batch) 的张量，分别对应水平极化(H)和垂直极化(V)的纯土壤发射亮温 (单位: K)
        """
        device = t_soi_k.device
        dtype = t_soi_k.dtype
        B, N = t_soi_k.shape

        # 1. 基础物理参数准备（挂在 Tensor 上）
        f = fghz * torch.tensor(1e9, device=device, dtype=dtype)
        lam = torch.as_tensor(self.C, device=device, dtype=dtype) / f
        k = torch.tensor(2.0, device=device, dtype=dtype) * torch.as_tensor(self.pi, device=device, dtype=dtype) / lam
        kcm = k / torch.tensor(100.0, device=device, dtype=dtype)

        # 2. 计算各层介电常数
        wf_clay_2d = wf_clay_surf.unsqueeze(1).expand(-1, N)
        f_2d = f.unsqueeze(1).expand(-1, N) if f.dim() > 0 else f
        eps_soil_profile_liq = self.diel_soil_M09(liq_soi_profile, t_soi_k, wf_clay_2d, f_2d)
        
        eps_f = torch.complex(torch.full_like(liq_soi_profile, 5.0), torch.full_like(liq_soi_profile, 0.5))
        eps_soil = eps_soil_profile_liq * (torch.tensor(1.0, device=device, dtype=dtype) - ffrz_profile) + eps_f * ffrz_profile

        eps_r = eps_soil.real
        eps_i = torch.abs(eps_soil.imag)

        # ================== 替换这里的代码 ==================
        # 3. 计算表层菲涅尔反射率及粗糙度修正
        eps_soil_surf = eps_soil[:, 0]
        
        # 修复：计算表层时直接使用 1D 张量，避免广播出 (15000, 15000) 的错位矩阵
        sin2_theta_surf = torch.sin(theta)**2
        cos_theta_surf = torch.cos(theta)

        g_val = torch.sqrt(eps_soil_surf - sin2_theta_surf)
        r_s_h = torch.abs((cos_theta_surf - g_val)/(cos_theta_surf + g_val))**2
        r_s_v = torch.abs((cos_theta_surf*eps_soil_surf - g_val)/(cos_theta_surf*eps_soil_surf + g_val))**2
        r_s = torch.stack([r_s_h, r_s_v])

        # 获取粗糙地表反射率
        r_r = self.rough_reflectivity(is_desert, patchclass, r_s, theta, fghz, kcm, self.rgh_surf)
        r0_h, r0_v = r_r[0], r_r[1]
        
        if r0_h.dim() == 1:
            r0_h = r0_h.unsqueeze(1)
            r0_v = r0_v.unsqueeze(1)

        # 新增：为第 4 步的多层矩阵运算单独准备 2D 的角度张量
        if theta.dim() == 0:
            sin2_theta = torch.sin(theta)**2
        else:
            sin2_theta = torch.sin(theta).view(-1, 1)**2

        # 4. 计算局部折射角余弦 mu_j
        mod_eps = torch.sqrt((eps_r - sin2_theta)**2 + eps_i**2)
        mu_j = torch.sqrt((mod_eps + eps_r - sin2_theta) / (mod_eps + eps_r + sin2_theta))

        # ================== 修改开始 ==================
        # 🌟 修复：将 1D 的波数 k 强制扩展为 2D (15000, 1)，以匹配 10 层的 eps_soil
        k_2d = k.unsqueeze(1) if k.dim() == 1 else k
        
        # 5. 吸收系数 ka_j (使用 k_2d)
        ka_j = torch.tensor(2.0, device=device, dtype=dtype) * k_2d * torch.imag(torch.sqrt(eps_soil))

        # 6. 体散射系数 ks_j 与不对称因子 g_j (基于密集介质理论 DMT)
        wf_sand_2d = wf_sand.unsqueeze(1) if wf_sand.dim() == 1 else wf_sand
        porsl_2d = porsl.unsqueeze(1) if porsl.dim() == 1 else porsl

        r_j = (torch.tensor(0.01, device=device, dtype=dtype) + (torch.tensor(0.5, device=device, dtype=dtype) - torch.tensor(0.01, device=device, dtype=dtype)) * wf_sand_2d) * torch.tensor(1e-3, device=device, dtype=dtype)
        f_j = (torch.tensor(1.0, device=device, dtype=dtype) - porsl_2d) * wf_sand_2d

        eps_grain = torch.complex(torch.tensor(4.0, device=device, dtype=dtype), torch.tensor(0.05, device=device, dtype=dtype))
        y_R = ((eps_grain - torch.tensor(1.0, device=device, dtype=dtype)) / (eps_grain + torch.tensor(2.0, device=device, dtype=dtype))).real

        term1 = (torch.tensor(1.0, device=device, dtype=dtype) - f_j)**4 / (torch.tensor(1.0, device=device, dtype=dtype) + torch.tensor(2.0, device=device, dtype=dtype) * f_j)**2
        term2 = (torch.tensor(1.0, device=device, dtype=dtype) - f_j * y_R)**1.5 * (torch.tensor(1.0, device=device, dtype=dtype) + torch.tensor(2.0, device=device, dtype=dtype) * y_R)**0.5
        
        # 🌟 修复：这里的 k 也必须使用 k_2d，否则 (15000,) * (15000, 1) 会错误广播成 (15000, 15000) 的内存炸弹
        Qs_j = (torch.tensor(8.0 / 3.0, device=device, dtype=dtype)) * (k_2d * r_j)**4 * (y_R**2) * (term1 / term2)

        ks_j = (torch.tensor(3.0, device=device, dtype=dtype) * f_j) / (torch.tensor(4.0, device=device, dtype=dtype) * r_j) * Qs_j
        g_j = torch.tensor(0.23, device=device, dtype=dtype) * (k_2d * r_j)**2
        # ================== 修改结束 ==================

        # 【调试 Print】：在进行易错的加减运算前，如果发现维度依旧不对称，立刻打印日志
        if ks_j.dim() == 1 and ka_j.dim() == 2:
            print(f"\n[Debug] 维度不匹配警告! ks_j shape: {ks_j.shape}, ka_j shape: {ka_j.shape}")

        # 7. 二流 RTE 参数
        eps_stable = torch.tensor(1e-12, device=device, dtype=dtype)
        omega_j = ks_j / (ks_j + ka_j + eps_stable)
        a_j = torch.sqrt((torch.tensor(1.0, device=device, dtype=dtype) - omega_j) / (torch.tensor(1.0, device=device, dtype=dtype) - omega_j * g_j + eps_stable))
        L_j = a_j / (ka_j + eps_stable)
        kappa_j = torch.tensor(1.0, device=device, dtype=dtype) / (mu_j * L_j + eps_stable)
        tau_j = torch.exp(-kappa_j * self.dz_soi)
        r_inf_j = (torch.tensor(1.0, device=device, dtype=dtype) - a_j) / (torch.tensor(1.0, device=device, dtype=dtype) + a_j + eps_stable)

        # 8. 计算层间物理折射项
        p_j = torch.sqrt(eps_soil - sin2_theta)
        
        # 预分配块三对角矩阵的求解函数
        def solve_matrix(R0, R_int):
            M = torch.zeros((B, 2*N, 2*N), dtype=dtype, device=device)
            C = torch.zeros((B, 2*N, 1), dtype=dtype, device=device)

            R0_sq = R0.squeeze()
            r_inf_j_0 = r_inf_j[:, 0]

            # 顶部边界条件
            M[:, 0, 0] = r_inf_j_0 - R0_sq
            M[:, 0, 1] = torch.tensor(1.0, device=device, dtype=dtype) - R0_sq * r_inf_j_0
            C[:, 0, 0] = -(torch.tensor(1.0, device=device, dtype=dtype) - R0_sq) * t_soi_k[:, 0]

            # 内部层间边界条件
            if N > 1:
                # 核心修正位置：直接使用整个 R_int 代表 N-1 个内部界面的反射率
                # 这样 R_curr, tau, r_inf 的第二维全为 9 (即 N-1)，彻底消除尺寸不匹配
                R_curr = R_int 
                tau = tau_j[:, :-1]
                tau_inv = torch.tensor(1.0, device=device, dtype=dtype) / (tau + eps_stable)
                r_inf = r_inf_j[:, :-1]
                r_inf_p1 = r_inf_j[:, 1:]
                delta_t = t_soi_k[:, 1:] - t_soi_k[:, :-1]

                # 构建全局并发网格索引
                j_idx = torch.arange(N - 1, device=device)
                row1 = 2 * j_idx + 1
                row2 = 2 * j_idx + 2
                col0 = 2 * j_idx
                col1 = 2 * j_idx + 1
                col2 = 2 * j_idx + 2
                col3 = 2 * j_idx + 3

                # 对常数使用统一的 Tensor 运算
                one_t = torch.tensor(1.0, device=device, dtype=dtype)

                # Eq 1: 向上辐射边界匹配
                M[:, row1, col0] = (one_t - R_curr * r_inf) * tau_inv
                M[:, row1, col1] = (r_inf - R_curr) * tau
                M[:, row1, col2] = -(one_t - R_curr)
                M[:, row1, col3] = -(one_t - R_curr) * r_inf_p1
                C[:, row1, 0] = (one_t - R_curr) * delta_t

                # Eq 2: 向下辐射边界匹配
                M[:, row2, col0] = -(one_t - R_curr) * r_inf * tau_inv
                M[:, row2, col1] = -(one_t - R_curr) * tau
                M[:, row2, col2] = r_inf_p1 - R_curr
                M[:, row2, col3] = one_t - R_curr * r_inf_p1
                C[:, row2, 0] = -(one_t - R_curr) * delta_t

            # 底部半无限大边界条件
            M[:, 2*N - 1, 2*N - 2] = torch.tensor(1.0, device=device, dtype=dtype)
            C[:, 2*N - 1, 0] = torch.tensor(0.0, device=device, dtype=dtype)

            # 求解 2Nx2N 线性方程组得到系数 A 和 B
            X = torch.linalg.solve(M, C)
            A1 = X[:, 0, 0]
            B1 = X[:, 1, 0]

            # 提取最终纯土壤发射亮温
            tb = (torch.tensor(1.0, device=device, dtype=dtype) - R0_sq) * (A1 + B1 * r_inf_j_0 + t_soi_k[:, 0])
            return tb

        # 9. 计算水平极化 (H) 的层间 Fresnel 反射率及解算亮温
        # p_j 有 10 层，层间计算切片后刚好产生 9 个界面值
        R_int_h = torch.abs((p_j[:, :-1] - p_j[:, 1:]) / (p_j[:, :-1] + p_j[:, 1:]))**2
        tb_soil_h = solve_matrix(r0_h, R_int_h)

        # 10. 计算垂直极化 (V) 的层间 Fresnel 反射率及解算亮温
        eps_p_j = eps_soil[:, 1:] * p_j[:, :-1]
        eps_p_jp1 = eps_soil[:, :-1] * p_j[:, 1:]
        R_int_v = torch.abs((eps_p_j - eps_p_jp1) / (eps_p_j + eps_p_jp1))**2
        tb_soil_v = solve_matrix(r0_v, R_int_v)

        return torch.stack([tb_soil_h, tb_soil_v])

    def calculate_tb_soil_LandEM(self, t_skin, t_soil, liq_surf, wf_clay_surf, ffrz_surf,
                                 theta, fghz, is_desert, patchclass, 
                                 C_2=0.0479924, emiss_default_h=0.25, emiss_default_v=0.30):
        """
        基于 ARMS LandEM 方案简化的裸土（Open Lands）微波亮度温度计算。
        
        针对裸土等开阔区域的应用场景，函数排除了植被、雪层的影响（即单次散射反照率 ssalb=0, 光学厚度 tau=0）。
        ARMS 原版的双流近似积分方程在此物理边界下发生了解析退化，可省略复杂的多次散射与衰减求解过程，
        退化后的有效发射率仅受土壤表面粗糙反射率和表层-深层土壤温度的普朗克（Rayleigh-Jeans）耦合因子的控制。
        接口已与 CMEM 方案对齐，直接输入土壤物理信息，内部自动完成介电常数及粗糙度反射率的求解。

        参数说明:
        t_skin           : [tensor] 地表皮肤物理温度，形状 (batch,)，单位: K
        t_soil           : [tensor] 深层有效土壤物理温度，形状 (batch,)，单位: K
        liq_surf         : [tensor] 表面液态水体积含水量，形状 (batch,)，单位: m3/m3
        wf_clay_surf     : [tensor] 表面粘土质量百分含量，形状 (batch,)，单位: %
        ffrz_surf        : [tensor] 表面冻结比例，形状 (batch,)，范围: 0~1
        theta            : [tensor/float] 观测天顶角，单位: rad
        fghz             : [tensor/float] 微波频率，单位: GHz
        is_desert        : [tensor] 是否为荒漠的布尔掩码，形状 (batch,)
        patchclass       : [tensor] 地表覆盖类型索引，形状 (batch,)
        C_2              : [float] 普朗克-玻尔兹曼常数项 hc/kB，常数固定为 0.0479924
        emiss_default_h  : [float] 水平极化下限截断保护值 (默认 0.25)
        emiss_default_v  : [float] 垂直极化下限截断保护值 (默认 0.30)

        返回:
        tb_soil          : [tensor] 形状为 (2, batch) 的张量，分别对应水平极化(H)和垂直极化(V)的土壤向上发射亮温 (单位: K)
        """
        # 1. 物理温度范围校验与回退机制
        invalid_soil = (t_soil <= 100.0) | (t_soil >= 350.0)
        valid_skin = (t_skin >= 100.0) & (t_skin <= 350.0)
        t_soil_eff = torch.where(invalid_soil & valid_skin, t_skin, t_soil)

        # 2. 计算表面介电常数
        f = fghz * 1e9
        kcm = (2 * self.pi / (self.C / f)) / 100.0
        
        eps_f = torch.complex(torch.full_like(liq_surf, 5.0), torch.full_like(liq_surf, 0.5))
        eps_soil_surf_liq = self.diel_soil_M09(liq_surf, t_skin, wf_clay_surf, f)
        eps_soil_surf = eps_soil_surf_liq * (1.0 - ffrz_surf) + eps_f * ffrz_surf

        # 3. 计算地表粗糙度修正后的反射率 (对应原方程中的 r23_h, r23_v)
        g_val = torch.sqrt(eps_soil_surf - torch.sin(theta)**2)
        r_s_h = torch.abs((torch.cos(theta) - g_val)/(torch.cos(theta) + g_val))**2
        r_s_v = torch.abs((torch.cos(theta)*eps_soil_surf - g_val)/(torch.cos(theta)*eps_soil_surf + g_val))**2
        r_s = torch.stack([r_s_h, r_s_v])

        r_r = self.rough_reflectivity(is_desert, patchclass, r_s, theta, fghz, kcm, self.rgh_surf)
        r23_h, r23_v = r_r[0], r_r[1]

        # 4. 计算地表温度与有效土壤温度的耦合因子
        term_skin = torch.exp(C_2 * fghz / t_skin) - 1.0
        term_soil = torch.exp(C_2 * fghz / t_soil_eff) - 1.0
        gsect0 = term_skin / term_soil

        # 5. 求解等效裸土发射率 (无植被衰减，方程退化简化的最终形式)
        esh = (1.0 - r23_h) * gsect0
        esv = (1.0 - r23_v) * gsect0

        # 6. 物理边界截断控制
        esh = torch.clamp(esh, min=emiss_default_h, max=1.0)
        esv = torch.clamp(esv, min=emiss_default_v, max=1.0)

        # 7. 计算向上发射亮温
        tb_soil_h = t_skin * esh
        tb_soil_v = t_skin * esv

        return torch.stack([tb_soil_h, tb_soil_v])
    
    # =========================================================================
    # 以下为内部物理模型函数
    # =========================================================================

    def rough_reflectivity(self, is_desert, patchclass, r_s, theta, fghz, kcm, rgh_surf):
        """
        计算地表粗糙度影响下的反射率
        完全参照 Fortran 版本中的查找表与分支逻辑实现
        """
        p_class = patchclass.long()  
        
        # 计算极化混合参数 Q
        Q = torch.where(fghz < 2.0, 
                        torch.zeros_like(fghz),                   
                        0.35 * (1.0 - torch.exp(-0.6 * (rgh_surf**2) * fghz)))  

        # 根据 DEF_DA_RTM_rough 选项配置粗糙度参数 hr, nrh, nrv
        if self.def_da_rtm_rough == 0:
            hr = (2.0 * kcm * rgh_surf)**2.0
            nrh = torch.zeros_like(r_s[0])
            nrv = torch.zeros_like(r_s[0])
        elif self.def_da_rtm_rough == 1:
            hr = self.hr_SMOS[p_class]
            nrh = torch.full_like(r_s[0], 2.0)
            nrv = torch.zeros_like(r_s[0])
        elif self.def_da_rtm_rough == 2:
            hr = self.hr_SMAP[p_class]
            nrh = torch.full_like(r_s[0], 2.0)
            nrv = torch.full_like(r_s[0], 2.0)
        elif self.def_da_rtm_rough == 3:
            hr = self.hr_P16[p_class]
            nrh = torch.full_like(r_s[0], -1.0)
            nrv = torch.full_like(r_s[0], -1.0)
        else:
            hr = torch.zeros_like(r_s[0])
            nrh = torch.zeros_like(r_s[0])
            nrv = torch.zeros_like(r_s[0])
            
        # 计算包含粗糙度修正后的水平与垂直极化反射率
        r_r_h = (Q * r_s[1] + (1.0 - Q) * r_s[0]) * torch.exp(-hr * (torch.cos(theta)**nrh))
        r_r_v = (Q * r_s[0] + (1.0 - Q) * r_s[1]) * torch.exp(-hr * (torch.cos(theta)**nrv))
        
        r_r_rough = torch.stack([r_r_h, r_r_v])
        
        # 荒漠地区直接使用平滑反射率
        is_desert_2d = is_desert.unsqueeze(0).expand(2, -1)
        r_r_final = torch.where(is_desert_2d, r_s, r_r_rough)

        return r_r_final

    def eff_soil_temp_Wilheit(self, t_soi_k, liq_soi_profile, wf_clay_surf, f, ffrz_profile, theta, lamcm):
        """Wilheit 1978 相干模型，内部计算介电常数廓线"""
        # 计算各层介电常数 (冻土混合)
        eps_f = torch.complex(torch.tensor(5.0, device=t_soi_k.device), torch.tensor(0.5, device=t_soi_k.device))
        nlay = liq_soi_profile.shape[1]
        wf_clay_2d = wf_clay_surf.unsqueeze(1).expand(-1, nlay)
        f_2d = f.unsqueeze(1).expand(-1, nlay) if f.dim()>0 else f
        eps_liq = self.diel_soil_M09(liq_soi_profile, t_soi_k, wf_clay_2d, f_2d)
        eps = eps_liq * (1.0 - ffrz_profile) + eps_f * ffrz_profile

        batch_size = eps.shape[0]
        nlay = eps.shape[1]
        N = nlay + 1
        device = eps.device

        dtype_c = torch.complex64 
        dtype_f = torch.float32

        CN = torch.zeros((batch_size, N), dtype=dtype_c, device=device)
        CN[:, 0] = 1.0 + 0.0j
        CN[:, 1:] = torch.sqrt(eps)

        DEL = torch.zeros((batch_size, N), dtype=dtype_f, device=device)
        DEL[:, 1:] = self.dz_soi * 100.0   # 使用内置 dz_soi，单位 cm

        S = torch.sin(theta)
        CP = torch.ones((batch_size, N), dtype=dtype_c, device=device)
        
        NMAX_idx = torch.full((batch_size,), N-1, dtype=torch.long, device=device)
        active_mask = torch.ones(batch_size, dtype=torch.bool, device=device)

        for i in range(1, N-1):
            CS = CN[:, 0] * S / CN[:, i]
            CC = torch.sqrt(1.0 + 0.0j - CS*CS)
            ARG = DEL[:, i] * 2.0 * self.pi / lamcm
            CARG = 2.0 * ARG * CN[:, i] * CC * 1.0j
            
            CP_next = torch.exp(CARG) * CP[:, i-1]
            CP[:, i] = torch.where(active_mask, CP_next, CP[:, i-1])
            
            below_thresh = torch.abs(CP[:, i]) < 0.0001
            crossed_now = active_mask & below_thresh
            
            NMAX_idx = torch.where(crossed_now, torch.tensor(i + 1, device=device), NMAX_idx)
            active_mask = active_mask & (~below_thresh)

        CEP_h = torch.zeros((batch_size, N), dtype=dtype_c, device=device)
        CEM_h = torch.zeros((batch_size, N), dtype=dtype_c, device=device)
        CEP_v = torch.zeros((batch_size, N), dtype=dtype_c, device=device)
        CEM_v = torch.zeros((batch_size, N), dtype=dtype_c, device=device)

        batch_indices = torch.arange(batch_size, device=device)
        CEP_h[batch_indices, NMAX_idx] = 1.0 + 0.0j
        CEP_v[batch_indices, NMAX_idx] = 1.0 + 0.0j

        for j in range(N-2, -1, -1):
            active_j = j < NMAX_idx  
            safe_cp = torch.where(active_j, CP[:, j], torch.tensor(1.0, dtype=dtype_c, device=device))

            CSJ = CN[:, 0] * S / CN[:, j]
            CCJ = torch.sqrt(1.0 + 0.0j - CSJ*CSJ)
            CSJP1 = CN[:, 0] * S / CN[:, j+1]
            CCJP1 = torch.sqrt(1.0 + 0.0j - CSJP1*CSJP1)

            CA_h = 2.0 * CN[:, j] * CCJ / (CN[:, j]*CCJ + CN[:, j+1]*CCJP1)
            CB_h = (CN[:, j]*CCJ - CN[:, j+1]*CCJP1) / ((CN[:, j]*CCJ + CN[:, j+1]*CCJP1) * safe_cp)
            
            CEP_h_j = CEP_h[:, j+1]/CA_h + CB_h*CEM_h[:, j+1]/CA_h
            CEM_h_j = CEM_h[:, j+1] + (CEP_h[:, j+1] - CEP_h_j)*safe_cp

            CEP_h[:, j] = torch.where(active_j, CEP_h_j, CEP_h[:, j])
            CEM_h[:, j] = torch.where(active_j, CEM_h_j, CEM_h[:, j])

            CD_v = 2.0 * CN[:, j] * CCJ
            CA_v = CN[:, j]*CCJP1 + CN[:, j+1]*CCJ
            CB_v = CN[:, j]*CCJP1 - CN[:, j+1]*CCJ
            
            CEP_v_j = CA_v*CEP_v[:, j+1]/CD_v + CB_v*CEM_v[:, j+1]/(CD_v*safe_cp)
            CR_v = CN[:, j+1] / CN[:, j]
            CEM_v_j = CR_v*CEM_v[:, j+1] + (CEP_v_j - CEP_v[:, j+1]*CR_v)*safe_cp
            
            CEP_v[:, j] = torch.where(active_j, CEP_v_j, CEP_v[:, j])
            CEM_v[:, j] = torch.where(active_j, CEM_v_j, CEM_v[:, j])

        CX_h = CEP_h[:, 0]
        CX_v = CEP_v[:, 0]
        safe_CX_h = torch.where(torch.abs(CX_h) > 0, CX_h, torch.tensor(1.0, dtype=dtype_c, device=device))
        safe_CX_v = torch.where(torch.abs(CX_v) > 0, CX_v, torch.tensor(1.0, dtype=dtype_c, device=device))
        
        for j in range(N):
            CEP_h[:, j] = CEP_h[:, j] / safe_CX_h
            CEM_h[:, j] = CEM_h[:, j] / safe_CX_h
            CEP_v[:, j] = CEP_v[:, j] / safe_CX_v
            CEM_v[:, j] = CEM_v[:, j] / safe_CX_v

        fa_h = torch.zeros((batch_size, nlay), dtype=dtype_f, device=device)
        fa_v = torch.zeros((batch_size, nlay), dtype=dtype_f, device=device)
        cos_theta = torch.cos(theta)

        j_indices = torch.arange(N, device=device).unsqueeze(0)
        mask_nmax = j_indices >= NMAX_idx.unsqueeze(1)
        CP = torch.where(mask_nmax, torch.tensor(1e-15j, dtype=dtype_c, device=device), CP)

        for j in range(1, N):
            active_j = j <= NMAX_idx  
            CS = S / CN[:, j]
            CC = torch.sqrt(1.0 + 0.0j - CS*CS)
            
            R = torch.abs(CP[:, j])
            S_abs = torch.abs(CP[:, j-1])
            
            safe_R = torch.where(R > 1e-25, R, torch.tensor(1e-25, device=device))
            safe_S = torch.where(S_abs > 1e-25, S_abs, torch.tensor(1e-25, device=device))

            E2_h = (safe_S - safe_R) * torch.abs(CEP_h[:, j])**2 + (1.0/safe_R - 1.0/safe_S) * torch.abs(CEM_h[:, j])**2
            DP_h = E2_h * (CN[:, j]*CC).real / cos_theta
            CXP_h = CEP_h[:, j] * torch.conj(CEM_h[:, j])
            X_h = 2.0 * ((CN[:, j]*CC / cos_theta).imag) * ( (CXP_h * CP[:, j-1]/safe_S).imag - (CXP_h * CP[:, j]/safe_R).imag )
            
            fa_h_val = DP_h - X_h
            fa_h[:, j-1] = torch.where(active_j, fa_h_val.real, torch.zeros_like(fa_h_val.real))
            
            E2_v = (safe_S - safe_R) * torch.abs(CEP_v[:, j])**2 + (1.0/safe_R - 1.0/safe_S) * torch.abs(CEM_v[:, j])**2
            DP_v = E2_v * (CN[:, j]*CC).real / cos_theta
            CXP_v = CEP_v[:, j] * torch.conj(CEM_v[:, j])
            X_v = 2.0 * ((CN[:, j]*CC / cos_theta).imag) * ( (CXP_v * CP[:, j-1]/safe_S).imag - (CXP_v * CP[:, j]/safe_R).imag )
            
            fa_v_val = DP_v - X_v
            fa_v[:, j-1] = torch.where(active_j, fa_v_val.real, torch.zeros_like(fa_v_val.real))

        sum_fa_h = torch.sum(fa_h, dim=1)
        teff_h = torch.where(sum_fa_h == 0, t_soi_k[:, 0], torch.sum(fa_h * t_soi_k, dim=1) / sum_fa_h)
        
        sum_fa_v = torch.sum(fa_v, dim=1)
        teff_v = torch.where(sum_fa_v == 0, t_soi_k[:, 0], torch.sum(fa_v * t_soi_k, dim=1) / sum_fa_v)

        r_h = torch.abs(CEM_h[:, 0])**2 * CN[:, 0].real
        r_v = torch.abs(CEM_v[:, 0])**2 * CN[:, 0].real
        
        return r_h, r_v, teff_h, teff_v
    
    def eff_soil_temp_Wigneron2001(self, wc_surf, t_surf, t_deep):
        """Wigneron 2001 scheme"""
        w0 = 0.30
        bw = 0.30
        wc_safe = torch.clamp(wc_surf, min=1e-6)
        C = torch.clamp((wc_safe / w0)**bw, min=0.001)
        teff = t_deep + (t_surf - t_deep) * C
        return teff
    
    def eff_soil_temp_Holmes2006(self, eps_surf, t_surf, t_deep):
        """Holmes 2006 scheme"""
        eps_r = eps_surf.real
        eps_i = torch.abs(eps_surf.imag)
        b_param = 0.05      
        eps0_param = 0.08
        eps_ratio = eps_i / eps_r
        C = (eps_ratio / eps0_param) ** b_param
        C = torch.clamp(C, min=0.001)
        teff = t_deep + (t_surf - t_deep) * C
        return teff

    def desert(self, t_soil, r_r, eps, kr):
        """处理荒漠情况的发射率与亮温"""
        f0 = 0.7
        y_r = (eps.real - 1) / (eps.real + 2)                                   
        y_i = 3 * eps.imag / (eps.real + 2)**2                                  
        w = ((1 - f0)**4 * kr**3 * y_r**2) / ((1 - f0)**4 * kr**3 * y_r**2 + 1.5*(1 + 2*f0)**2 * y_i) 
        g = 0.23 * kr**2                                                        
        a = torch.sqrt((1 - w) / (1 - w*g))                                     
        em = (1 - r_r) * (2*a / ((1 + a) - (1 - a)*r_r))                        
        return t_soil * em