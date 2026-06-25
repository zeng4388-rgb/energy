#!/usr/bin/env python3
"""FRB Catalog 预处理脚本

读取原始 FITS catalog，按筛选条件过滤，输出干净数据。
分析代码统一读取预处理后的数据，不直接接触原始 catalog。

用法:
    python preprocess_catalog.py
    python preprocess_catalog.py --input data/chimefrbcat2.fits --output data/filtered
"""

import numpy as np
import argparse
from pathlib import Path
from astropy.io import fits


def load_raw_catalog(fname):
    """读取原始 FITS catalog"""
    with fits.open(fname) as hdul:
        data = hdul[1].data
        col_names = data.columns.names
        print(f"原始数据: {len(data)} 个 burst")
        print(f"列名: {col_names}")
        return data, col_names


def filter_minimum(data):
    """最小筛选：只排除重复暴和无fluence的burst

    筛选标准:
    1. repeater_name 为空 → one-off FRB
    2. fluence 不为 NaN → 有有效通量测量
    """
    n_raw = len(data)

    # 排除重复暴
    rep = np.array(data['repeater_name'])
    is_oneoff = np.array([str(r).strip() == '' for r in rep])
    n_after_rep = np.sum(is_oneoff)

    # 排除无fluence
    fluence = np.array(data['fluence'], dtype=float)
    has_fluence = ~np.isnan(fluence)
    n_after_fluence = np.sum(is_oneoff & has_fluence)

    # 组合筛选
    mask = is_oneoff & has_fluence
    filtered = data[mask]

    print(f"\n筛选结果:")
    print(f"  原始:        {n_raw}")
    print(f"  排除重复暴:  {n_raw} → {n_after_rep} (排除 {n_raw - n_after_rep})")
    print(f"  排除无fluence: {n_after_rep} → {n_after_fluence} (排除 {n_after_rep - n_after_fluence})")
    print(f"  最终样本:    {n_after_fluence} 个 one-off FRB")

    return filtered


def extract_columns(data):
    """提取分析所需列，统一单位

    输出列:
    - name: FRB 名称
    - fluence: 通量×持续时间 [Jy·ms]
    - width: 脉冲宽度 [ms]（从秒转换）
    - dm_obs: 观测色散量 [pc·cm⁻³]
    - dm_exc_ne2001: DM - DM_MW(NE2001) [pc·cm⁻³]
    - dm_exc_ymw16: DM - DM_MW(YMW16) [pc·cm⁻³]
    """
    n = len(data)
    result = {}

    # 名称
    result['name'] = np.array(data['tns_name'])

    # fluence [Jy·ms]
    result['fluence'] = np.array(data['fluence'], dtype=float)

    # width: CHIME bc_width 单位是秒，转为毫秒
    width = np.array(data['bc_width'], dtype=float)
    med_w = np.nanmedian(width)
    if np.isfinite(med_w) and med_w < 1.0:
        width = width * 1000.0
    result['width'] = width

    # DM
    result['dm_obs'] = np.array(data['dm_fitb'], dtype=float)
    result['dm_exc_ne2001'] = np.array(data['dm_exc_ne2001'], dtype=float)
    result['dm_exc_ymw16'] = np.array(data['dm_exc_ymw16'], dtype=float)

    return result


def save_fits(result, fname):
    """保存为 FITS 格式"""
    cols = []
    for key in result:
        if result[key].dtype.kind in ('U', 'S'):
            # 字符串列：使用最大长度避免截断
            max_len = max(len(s) for s in result[key])
            col = fits.Column(name=key, format=f'{max_len}A', array=result[key])
        else:
            col = fits.Column(name=key, format='D', array=result[key])
        cols.append(col)

    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.header['CATALOG'] = 'CHIME/FRB Catalog 2 (filtered)'
    hdu.header['NAXIS2'] = len(result['fluence'])
    hdu.header['WUNIT'] = 'ms'

    hdul = fits.HDUList([fits.PrimaryHDU(), hdu])
    hdul.writeto(fname, overwrite=True)
    print(f"  FITS: {fname}")


def save_txt(result, fname):
    """保存为 TXT 格式（空格分隔，首行为列名）"""
    keys = list(result.keys())
    header = ' '.join(keys)

    # 分离字符串列和数值列
    str_keys = [k for k in keys if result[k].dtype.kind in ('U', 'S', 'O')]
    num_keys = [k for k in keys if result[k].dtype.kind not in ('U', 'S', 'O')]

    # 用混合格式写入
    with open(fname, 'w') as f:
        f.write(f'# {header}\n')
        for i in range(len(result[keys[0]])):
            parts = []
            for k in keys:
                val = result[k][i]
                if k in str_keys:
                    parts.append(str(val))
                else:
                    parts.append(f'{float(val):.6f}')
            f.write(' '.join(parts) + '\n')
    print(f"  TXT:  {fname}")


def print_summary(result):
    """打印数据摘要"""
    print(f"\n数据摘要:")
    for key in result:
        vals = result[key]
        if vals.dtype.kind in ('U', 'S', 'O'):
            print(f"  {key:>15}: {len(vals)} 条")
        else:
            print(f"  {key:>15}: min={np.nanmin(vals):.4f}  max={np.nanmax(vals):.4f}  median={np.nanmedian(vals):.4f}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='FRB Catalog 预处理')
    parser.add_argument('-i', '--input', default='data/chimefrbcat2.fits',
                        help='输入 FITS 文件路径')
    parser.add_argument('-o', '--output', default='data/filtered',
                        help='输出文件路径（不含扩展名）')
    args = parser.parse_args()

    print("=" * 50)
    print("FRB Catalog 预处理")
    print("=" * 50)

    # 1. 加载原始数据
    print(f"\n[1/4] 加载原始数据: {args.input}")
    data, col_names = load_raw_catalog(args.input)

    # 2. 筛选
    print(f"\n[2/4] 最小筛选（排除重复暴 + 无fluence）")
    filtered = filter_minimum(data)

    # 3. 提取列并转换单位
    print(f"\n[3/4] 提取列并转换单位")
    result = extract_columns(filtered)
    print_summary(result)

    # 4. 保存
    print(f"\n[4/4] 保存预处理数据")
    fits_out = args.output + '.fits'
    txt_out = args.output + '.txt'
    Path(fits_out).parent.mkdir(parents=True, exist_ok=True)
    save_fits(result, fits_out)
    save_txt(result, txt_out)

    print(f"\n完成！分析代码使用预处理后的数据:")
    print(f'  config.json: "catalog_path": "{fits_out}"')