#!/home/albatross/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
- TRC  : Total Read Count
- FF_snp: SNP-based Fetal Fraction
- XY_sample_df: Dataframe for unique_file
- ACA_ref_df / Y_ref_df: reference for ACA Data Generation(.unique file)
- C_yRC' (beta0 + beta1 * TRC), C_yRC (= C_yRC' + ε, ε~N(0, σ^2))
- Equation (1) ACA: T_iRC = C_iRC + (C_iRC * FF / 2)
- Equation (3) XO : C_XRC = C_XYRC - C_YRC + C_yRC
- Equation (4) XYY: C_XYYRC = C_XYRC + C_YRC - C_yRC
"""
import sys
import numpy as np
import pandas as pd


# ───────────────────────────────────────────────────────────
# Define bins (300kb) 
# ───────────────────────────────────────────────────────────
BIN_SIZE = 300_000
# chr location with terminal
CHR_END = {
    1: 249_300_000,  2: 243_300_000,  3: 198_300_000,  4: 191_400_000,
    5: 181_200_000,  6: 171_300_000,  7: 159_300_000,  8: 146_400_000,
    9: 141_300_000, 10: 135_600_000, 11: 135_300_000, 12: 134_100_000,
    13: 115_200_000, 14: 107_400_000, 15: 102_600_000, 16:  90_600_000,
    17:  81_300_000, 18:  78_300_000, 19:  59_400_000, 20:  63_300_000,
    21:  48_300_000, 22:  51_600_000, 23: 155_400_000, 24:  59_400_000,
}
# chr bin boundary lists
cnv_chr_range = {i: list(np.arange(0, CHR_END[i] + 1, BIN_SIZE)) for i in range(1, 25)}


# ───────────────────────────────────────────────────────────
# 유틸: 300kb Region count & generationg chr_* column name
# ───────────────────────────────────────────────────────────
def count_bins(df: pd.DataFrame) -> list:
    counts = []
    for i in range(1, 25):
        temp = df[df['Chr'] == f'chr{i}']
        edges = cnv_chr_range[i]
        for j in range(len(edges) - 1):
            left, right = edges[j], edges[j + 1]
            # 원 코드와 동일하게 좌/우 경계 모두 포함
            cnt = len(temp[(left <= temp['Pos']) & (temp['Pos'] <= right)])
            counts.append(cnt)
    return counts


def chr_colnames() -> list:
    names = []
    for i in range(1, 25):
        for j in range(len(cnv_chr_range[i]) - 1):
            names.append(f"chr{i}_{j+1}")
    return names


# # ───────────────────────────────────────────────────────────
# # Normal (음성) 샘플 CSV 생성
# # ───────────────────────────────────────────────────────────
# def make_normal():
#     count_read_pos = count_bins(XY_sample_df)
#     chr_cols = chr_colnames()

#     final_col_names = ['sample_id', 'result', 'GC', 'UR', 'snp_FF', 'del_range'] + chr_cols
#     final_df = pd.DataFrame(
#         np.array([file_names, 'M', Sample_GC, TRC, FF_snp, 0] + count_read_pos)
#     ).T
#     final_df.columns = final_col_names
#     final_df.to_csv(output_dir + file_names + '.csv', index=False)


# ───────────────────────────────────────────────────────────
# ACA(Trisomy) 합성 CSV 생성 — 식 (1): T_iRC = C_iRC + (C_iRC * FF / 2)
# ───────────────────────────────────────────────────────────
def make_ACA():
    ACA_ref_df = pd.read_table(working_dir + 'ACA_ref.unique', header=None)
    ACA_ref_df = ACA_ref_df[[2, 3, 9]]
    ACA_ref_df.columns = ['Chr', 'Pos', 'Seq']

    T_iRC_list = []
    S_ACA_df = XY_sample_df.copy()

    for i in range(1, 23):  # autosomes
        C_iRC = len(XY_sample_df[XY_sample_df['Chr'] == f'chr{i}'])
        delta_T_iRC = int(C_iRC * (FF_snp * 0.01) * 0.5)  # C_iRC * FF/2
        add_df = ACA_ref_df[ACA_ref_df['Chr'] == f'chr{i}'].sample(delta_T_iRC)
        S_ACA_df = pd.concat([S_ACA_df, add_df], axis=0, ignore_index=True)

        # 원 코드의 ACA_UR_list 계산 관습 유지(ATRC + 추가량)
        T_iRC = ATRC + delta_T_iRC
        T_iRC_list.append(T_iRC)

    count_read_pos = count_bins(S_ACA_df)
    chr_cols = chr_colnames()
    ACA_UR_col = [f"T{i}_UR" for i in range(1, 23)]

    final_col_names = ['sample_id', 'result', 'GC', 'UR', 'snp_FF'] + ACA_UR_col + chr_cols
    final_df = pd.DataFrame(
        np.array([file_names, 'ACA', Sample_GC, TRC, FF_snp] + T_iRC_list + count_read_pos)
    ).T
    final_df.columns = final_col_names
    final_df.to_csv(output_dir + file_names + '_ACA.csv', index=False)


# # ───────────────────────────────────────────────────────────
# # SCA(XO, XYY) 합성 CSV 생성
# # XO  : 식 (3), (3.1), (3.2) — 회귀 + 정규잔차로 C_yRC 모델링
# # XYY : 식 (4) 아이디어. 아래 구현은 '추가만' 근사(원 코드 유지).
# # ───────────────────────────────────────────────────────────
# def make_SCA():
#     Y_ref_df = pd.read_table(working_dir + 'Y_ref.unique', header=None)
#     Y_ref_df = Y_ref_df[[2, 3, 9]]
#     Y_ref_df.columns = ['Chr', 'Pos', 'Seq']
#     Y_ref_df['Chr'] = Y_ref_df['Chr'].replace(['chrX', 'chrY'], ['chr23', 'chr24'])

#     # 회귀 계수(식 3.2) & 잔차 표준편차(식 3.1)
#     beta0 = 8.73854394
#     beta1 = 3.76076E-05
#     sigma = 14

#     for disease_name in ['XO', 'XYY']:
#         if disease_name == 'XO':
#             # C_yRC' = β0 + β1·TRC (식 3.2)
#             C_yRC_pred = beta0 + beta1 * TRC
#             # C_yRC = C_yRC' + ε, ε~N(0, σ^2) (식 3.1)
#             C_yRC = int(np.random.normal(C_yRC_pred, sigma))

#             nonY_df = XY_sample_df[XY_sample_df['Chr'] != 'chr24']
#             Y_reads_XY_df = XY_sample_df[XY_sample_df['Chr'] == 'chr24']

#             # XO: 정상 Y 제거, 미스얼라인 C_yRC만 남김(원 코드 로직 유지)
#             S_SCA_df = pd.concat([nonY_df, Y_reads_XY_df.sample(C_yRC)], axis=0, ignore_index=True)

#         elif disease_name == 'XYY':
#             # XYY: 식 (4) 개념(C_XYYRC = C_XYRC + C_YRC - C_yRC)
#             # 현 구현은 '추가만' 방식(근사). 필요시 '제거+추가'로 바꿀 수 있음.
#             Y_reads_XY_df = XY_sample_df[XY_sample_df['Chr'] == 'chr24']

#             C_yRC_pred = beta0 + beta1 * TRC
#             C_yRC = int(np.random.normal(C_yRC_pred, sigma))

#             # 추가할 정상 Y 추정량(= C_XYRC - C_yRC)
#             C_YRC_est = max(0, len(Y_reads_XY_df) - C_yRC)

#             S_SCA_df = pd.concat([XY_sample_df, Y_ref_df.sample(C_YRC_est)], axis=0, ignore_index=True)

#         # 공통: 300kb bin 카운트 및 CSV 저장
#         count_read_pos = count_bins(S_SCA_df)
#         chr_cols = chr_colnames()
#         final_col_names = ['sample_id', 'result', 'GC', 'UR', 'snp_FF', 'del_range'] + chr_cols
#         final_df = pd.DataFrame(
#             np.array([file_names, disease_name, Sample_GC, TRC, FF_snp, 0] + count_read_pos)
#         ).T
#         final_df.columns = final_col_names
#         final_df.to_csv(output_dir + file_names + '_' + disease_name + '.csv', index=False)


# ───────────────────────────────────────────────────────────
# main
# ───────────────────────────────────────────────────────────
if __name__ == "__main__":
    files = sys.argv
    if len(files) < 2:
        print("Usage: script.py <Fastq_ID>")
        sys.exit(1)

    file_names = files[1]

    # Setting Diretory
    working_dir = '/BiO/snp_imp/jeong/test/'
    output_dir = '/BiO/snp_imp/jeong/test/'

    # Loading FF_snp
    FF_snp = float(file_names.split('_')[-2])

    # Load Reads info from .unique file 
    XY_sample_df = pd.read_table(working_dir + file_names + '.unique', header=None)
    XY_sample_df = XY_sample_df[[2, 3, 9]]
    XY_sample_df.columns = ['Chr', 'Pos', 'Seq']

    # Calculating GC Contents / Read Count
    XY_sample_df["GC"] = XY_sample_df["Seq"].apply(lambda x: x.count("G") + x.count("C"))
    XY_sample_df["Total"] = XY_sample_df["Seq"].apply(len)
    Sample_GC = XY_sample_df["GC"].sum() / XY_sample_df["Total"].sum()

    # chrX/Y → chr23/24
    XY_sample_df['Chr'] = XY_sample_df['Chr'].replace(['chrX', 'chrY'], ['chr23', 'chr24'])

    # TRC
    TRC = XY_sample_df.shape[0]
    # ATRC: Total Read Count of Autosomes
    ATRC = XY_sample_df[(XY_sample_df['Chr'] != 'chr23') & (XY_sample_df['Chr'] != 'chr24')].shape[0]

    # Fuction Calling
    make_ACA()
