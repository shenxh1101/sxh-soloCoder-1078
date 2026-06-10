#!/usr/bin/env python3
"""
测试新功能的脚本
1. 批量比对CSV/TSV导出
2. 保守性分析详细结果导出
3. DNA反向互补比对
4. GenBank缓存功能
"""

import subprocess
import sys
import os
import shutil

PYTHON = sys.executable
SCRIPT = "seq_align.py"
TEST_DIR = "test_data"
OUTPUT_DIR = "test_output"

def run_command(cmd, description):
    """运行命令并打印结果"""
    print(f"\n{'='*60}")
    print(f"测试: {description}")
    print(f"命令: {cmd}")
    print('='*60)
    
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
        if result.stdout:
            print("STDOUT:")
            print(result.stdout)
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
        print(f"返回码: {result.returncode}")
        return result.returncode == 0
    except Exception as e:
        print(f"执行出错: {e}")
        return False

def test_reverse_complement():
    """测试DNA反向互补比对功能"""
    print("\n" + "#"*60)
    print("# 测试1: DNA反向互补比对")
    print("#"*60)
    
    seq1 = os.path.join(TEST_DIR, "seq1.fasta")
    seq2 = os.path.join(TEST_DIR, "seq2_rc.fasta")
    
    # 不使用反向互补选项（应该得分较低）
    cmd1 = f"{PYTHON} {SCRIPT} -f {seq1} -s {seq2} -a global --no-color"
    result1 = run_command(cmd1, "不使用反向互补选项")
    
    # 使用反向互补选项（应该选择反向互补方向，得分较高）
    cmd2 = f"{PYTHON} {SCRIPT} -f {seq1} -s {seq2} -a global --reverse-complement --no-color"
    result2 = run_command(cmd2, "使用反向互补选项")
    
    return result1 and result2

def test_batch_export():
    """测试批量比对CSV/TSV导出"""
    print("\n" + "#"*60)
    print("# 测试2: 批量比对CSV/TSV导出")
    print("#"*60)
    
    ref = os.path.join(TEST_DIR, "batch", "ref.fasta")
    batch_dir = os.path.join(TEST_DIR, "batch")
    
    # 测试TXT格式导出
    txt_output = os.path.join(OUTPUT_DIR, "batch_results.txt")
    cmd1 = f"{PYTHON} {SCRIPT} --dir {batch_dir} --ref {ref} -a global -o {txt_output} --batch-format txt --no-color"
    result1 = run_command(cmd1, "批量比对TXT导出")
    
    # 测试CSV格式导出
    csv_output = os.path.join(OUTPUT_DIR, "batch_results.csv")
    cmd2 = f"{PYTHON} {SCRIPT} --dir {batch_dir} --ref {ref} -a global -o {csv_output} --batch-format csv --no-color"
    result2 = run_command(cmd2, "批量比对CSV导出")
    
    # 测试TSV格式导出
    tsv_output = os.path.join(OUTPUT_DIR, "batch_results.tsv")
    cmd3 = f"{PYTHON} {SCRIPT} --dir {batch_dir} --ref {ref} -a global -o {tsv_output} --batch-format tsv --no-color"
    result3 = run_command(cmd3, "批量比对TSV导出")
    
    # 验证文件存在
    for f in [txt_output, csv_output, tsv_output]:
        if os.path.exists(f):
            print(f"✓ 文件已生成: {f}")
            with open(f, 'r', encoding='utf-8') as fp:
                first_line = fp.readline().strip()
                print(f"  首行: {first_line}")
        else:
            print(f"✗ 文件未生成: {f}")
    
    return result1 and result2 and result3

def test_conservation_export():
    """测试保守性分析详细结果导出"""
    print("\n" + "#"*60)
    print("# 测试3: 保守性分析详细结果导出")
    print("#"*60)
    
    msa_files = [
        os.path.join(TEST_DIR, "msa", "seq_a.fasta"),
        os.path.join(TEST_DIR, "msa", "seq_b.fasta"),
        os.path.join(TEST_DIR, "msa", "seq_c.fasta"),
        os.path.join(TEST_DIR, "msa", "seq_d.fasta")
    ]
    
    # 测试CSV格式导出
    csv_output = os.path.join(OUTPUT_DIR, "conservation_details.csv")
    cmd1 = f"{PYTHON} {SCRIPT} --msa {' '.join(msa_files)} --conservation --conservation-export {csv_output} --conservation-format csv --no-color"
    result1 = run_command(cmd1, "保守性分析CSV导出")
    
    # 测试TSV格式导出
    tsv_output = os.path.join(OUTPUT_DIR, "conservation_details.tsv")
    cmd2 = f"{PYTHON} {SCRIPT} --msa {' '.join(msa_files)} --conservation --conservation-export {tsv_output} --conservation-format tsv --no-color"
    result2 = run_command(cmd2, "保守性分析TSV导出")
    
    # 测试TXT格式导出
    txt_output = os.path.join(OUTPUT_DIR, "conservation_details.txt")
    cmd3 = f"{PYTHON} {SCRIPT} --msa {' '.join(msa_files)} --conservation --conservation-export {txt_output} --conservation-format txt --no-color"
    result3 = run_command(cmd3, "保守性分析TXT导出")
    
    # 验证文件存在
    for f in [csv_output, tsv_output, txt_output]:
        if os.path.exists(f):
            print(f"✓ 文件已生成: {f}")
            with open(f, 'r', encoding='utf-8') as fp:
                lines = fp.readlines()
                print(f"  行数: {len(lines)}")
                if len(lines) > 1:
                    print(f"  表头: {lines[1].strip() if lines[0].startswith('#') else lines[0].strip()}")
        else:
            print(f"✗ 文件未生成: {f}")
    
    return result1 and result2 and result3

def test_cache_commands():
    """测试GenBank缓存管理命令"""
    print("\n" + "#"*60)
    print("# 测试4: GenBank缓存管理命令")
    print("#"*60)
    
    # 测试列出缓存
    cmd1 = f"{PYTHON} {SCRIPT} --list-cache"
    result1 = run_command(cmd1, "列出缓存")
    
    return result1

def test_help():
    """测试帮助信息是否包含新参数"""
    print("\n" + "#"*60)
    print("# 测试5: 帮助信息检查")
    print("#"*60)
    
    cmd = f"{PYTHON} {SCRIPT} -h"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    help_text = result.stdout
    new_params = [
        '--batch-format',
        '--conservation-export',
        '--conservation-format',
        '--reverse-complement',
        '--no-cache',
        '--clear-cache',
        '--list-cache'
    ]
    
    print("检查新参数是否在帮助中:")
    all_found = True
    for param in new_params:
        if param in help_text:
            print(f"  ✓ {param}")
        else:
            print(f"  ✗ {param} - 未找到!")
            all_found = False
    
    return all_found

def main():
    """主测试函数"""
    print("="*60)
    print("生物序列比对工具 - 新功能测试")
    print("="*60)
    
    # 创建输出目录
    if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)
    
    results = {}
    
    # 运行所有测试
    results['help'] = test_help()
    results['reverse_complement'] = test_reverse_complement()
    results['batch_export'] = test_batch_export()
    results['conservation_export'] = test_conservation_export()
    results['cache_commands'] = test_cache_commands()
    
    # 打印总结
    print("\n" + "="*60)
    print("测试总结")
    print("="*60)
    
    all_passed = True
    for test_name, passed in results.items():
        status = "✓ 通过" if passed else "✗ 失败"
        print(f"{test_name:30}: {status}")
        if not passed:
            all_passed = False
    
    print("\n" + "="*60)
    if all_passed:
        print("所有测试通过!")
    else:
        print("部分测试失败，请检查!")
    print("="*60)
    
    return 0 if all_passed else 1

if __name__ == '__main__':
    sys.exit(main())
