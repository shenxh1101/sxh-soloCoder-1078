#!/usr/bin/env python3
"""
生物序列比对可视化工具 - 主入口
支持DNA/蛋白质序列的全局和局部比对，彩色高亮输出，批量比对等功能
"""

import argparse
import sys
import os
from seq_align_tool import (
    read_fasta,
    read_fasta_directory,
    needleman_wunsch,
    smith_waterman,
    ScoringMatrix,
    color_alignment,
    print_alignment,
    export_alignment,
    batch_alignment,
    export_batch_results,
    conservation_analysis,
    print_conservation_result,
    export_conservation_details,
    fetch_genbank,
    reverse_complement,
    align_with_reverse_complement,
    AlignmentResult
)
from seq_align_tool.genbank import clear_cache, list_cache


def parse_args():
    parser = argparse.ArgumentParser(
        description='生物序列比对可视化工具 - DNA/Protein Sequence Alignment Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例用法:
  # 全局比对两个FASTA文件
  python seq_align.py -f seq1.fasta -s seq2.fasta -a global
  
  # 局部比对，自定义打分
  python seq_align.py -f seq1.fasta -s seq2.fasta -a local --match 2 --mismatch -1 --gap -2
  
  # 使用BLOSUM62矩阵进行蛋白质比对
  python seq_align.py -f prot1.fasta -s prot2.fasta -a global --matrix blosum62
  
  # 批量比对目录中的序列与参考序列
  python seq_align.py --dir ./sequences --ref reference.fasta -a global
  
  # 多序列保守性分析
  python seq_align.py --msa seq1.fasta seq2.fasta seq3.fasta --conservation
  
  # 从GenBank获取序列并比对
  python seq_align.py --genbank NM_001301717 --genbank2 NM_001301718 -a global
  
  # 导出结果为ALN格式
  python seq_align.py -f seq1.fasta -s seq2.fasta -a global -o result.aln --format aln
        '''
    )
    
    input_group = parser.add_argument_group('输入选项')
    input_group.add_argument('-f', '--file1', help='第一个FASTA序列文件')
    input_group.add_argument('-s', '--file2', help='第二个FASTA序列文件')
    input_group.add_argument('--dir', help='包含FASTA文件的目录（批量比对）')
    input_group.add_argument('--ref', help='批量比对时的参考序列文件')
    input_group.add_argument('--msa', nargs='+', help='多个FASTA文件用于保守性分析')
    input_group.add_argument('--genbank', help='GenBank accession号（序列1）')
    input_group.add_argument('--genbank2', help='GenBank accession号（序列2）')
    
    align_group = parser.add_argument_group('比对选项')
    align_group.add_argument('-a', '--algorithm', choices=['global', 'local'], default='global',
                           help='比对算法: global (Needleman-Wunsch) 或 local (Smith-Waterman)')
    align_group.add_argument('--match', type=int, default=1, help='匹配得分 (默认: 1)')
    align_group.add_argument('--mismatch', type=int, default=-1, help='错配得分 (默认: -1)')
    align_group.add_argument('--gap', type=int, default=-2, help='gap罚分 (默认: -2)')
    align_group.add_argument('--matrix', choices=['blosum45', 'blosum62', 'blosum80', 'pam30', 'pam70', 'pam250'],
                           help='蛋白质打分矩阵')
    align_group.add_argument('--matrix-file', help='自定义打分矩阵文件')
    align_group.add_argument('--seq-type', choices=['dna', 'protein', 'auto'], default='auto',
                           help='序列类型 (默认: auto自动检测)')
    align_group.add_argument('--reverse-complement', action='store_true',
                            help='DNA序列自动比较原序列和反向互补序列，选择最佳方向')
    
    output_group = parser.add_argument_group('输出选项')
    output_group.add_argument('-o', '--output', help='输出文件路径')
    output_group.add_argument('--format', choices=['txt', 'aln'], default='txt',
                             help='双序列比对输出格式 (默认: txt)')
    output_group.add_argument('--batch-format', choices=['txt', 'csv', 'tsv'], default='txt',
                             help='批量比对输出格式 (默认: txt)')
    output_group.add_argument('--no-color', action='store_true', help='禁用彩色输出')
    output_group.add_argument('--line-width', type=int, default=60, help='输出行宽 (默认: 60)')
    
    analysis_group = parser.add_argument_group('分析选项')
    analysis_group.add_argument('--conservation', action='store_true', help='进行保守性分析')
    analysis_group.add_argument('--plot-width', type=int, default=80, help='保守性曲线图宽度 (默认: 80)')
    analysis_group.add_argument('--conservation-export', help='保守性分析详细结果导出路径')
    analysis_group.add_argument('--conservation-format', choices=['txt', 'csv', 'tsv'], default='csv',
                               help='保守性分析详细结果格式 (默认: csv)')
    
    cache_group = parser.add_argument_group('缓存选项')
    cache_group.add_argument('--no-cache', action='store_true', help='禁用GenBank本地缓存')
    cache_group.add_argument('--clear-cache', action='store_true', help='清理GenBank本地缓存后退出')
    cache_group.add_argument('--list-cache', action='store_true', help='列出所有缓存的GenBank序列后退出')
    
    return parser.parse_args()


def detect_seq_type(seq):
    dna_chars = set('ATCGN')
    seq_upper = seq.upper().replace('-', '')
    if all(c in dna_chars for c in seq_upper):
        return 'dna'
    return 'protein'


def get_sequences(args):
    seq1, seq2, name1, name2 = None, None, None, None
    use_cache = not args.no_cache
    
    if args.genbank:
        print(f"正在从GenBank获取 {args.genbank}...", file=sys.stderr)
        name1, seq1 = fetch_genbank(args.genbank, use_cache=use_cache)
        if not seq1:
            print(f"错误: 无法获取GenBank序列 {args.genbank}", file=sys.stderr)
            sys.exit(1)
    
    if args.genbank2:
        print(f"正在从GenBank获取 {args.genbank2}...", file=sys.stderr)
        name2, seq2 = fetch_genbank(args.genbank2, use_cache=use_cache)
        if not seq2:
            print(f"错误: 无法获取GenBank序列 {args.genbank2}", file=sys.stderr)
            sys.exit(1)
    
    if args.file1:
        seqs = read_fasta(args.file1)
        if seqs:
            name1, seq1 = seqs[0]
    
    if args.file2:
        seqs = read_fasta(args.file2)
        if seqs:
            name2, seq2 = seqs[0]
    
    return seq1, seq2, name1, name2


def main():
    args = parse_args()
    
    if args.clear_cache:
        print("正在清理GenBank本地缓存...", file=sys.stderr)
        count = clear_cache()
        print(f"已清理 {count} 个缓存文件", file=sys.stderr)
        return
    
    if args.list_cache:
        cached = list_cache()
        if not cached:
            print("当前没有缓存的序列", file=sys.stderr)
        else:
            print(f"共有 {len(cached)} 个缓存的序列:")
            print(f"{'Accession':<20}{'数据库':<10}{'名称':<30}{'缓存时间':<20}")
            print("-" * 80)
            for item in cached:
                print(f"{item['accession']:<20}{item['db']:<10}{item['name']:<30}{item['cached_at']:<20}")
        return
    
    if args.msa and args.conservation:
        if len(args.msa) < 2:
            print("错误: 保守性分析需要至少2个序列", file=sys.stderr)
            sys.exit(1)
        
        sequences = []
        for fasta_file in args.msa:
            seqs = read_fasta(fasta_file)
            if seqs:
                sequences.extend(seqs)
        
        if len(sequences) < 2:
            print("错误: 未能读取到足够的序列", file=sys.stderr)
            sys.exit(1)
        
        print("\n" + "="*60)
        print("多序列保守性分析")
        print("="*60)
        
        result = conservation_analysis(sequences, plot_width=args.plot_width)
        print_conservation_result(result)
        
        if args.conservation_export:
            export_conservation_details(result, args.conservation_export, 
                                       fmt=args.conservation_format)
            print(f"\n保守性分析详细结果已导出到: {args.conservation_export}")
        return
    
    if args.dir and args.ref:
        ref_seqs = read_fasta(args.ref)
        if not ref_seqs:
            print(f"错误: 无法读取参考序列文件 {args.ref}", file=sys.stderr)
            sys.exit(1)
        
        ref_name, ref_seq = ref_seqs[0]
        query_seqs = read_fasta_directory(args.dir)
        
        if not query_seqs:
            print(f"错误: 目录 {args.dir} 中没有找到FASTA文件", file=sys.stderr)
            sys.exit(1)
        
        seq_type = args.seq_type
        if seq_type == 'auto':
            seq_type = detect_seq_type(ref_seq)
        
        scoring = ScoringMatrix(
            seq_type=seq_type,
            match=args.match,
            mismatch=args.mismatch,
            gap=args.gap,
            matrix_name=args.matrix,
            matrix_file=args.matrix_file
        )
        
        print("\n" + "="*80)
        print(f"批量比对 - 参考序列: {ref_name}")
        print(f"算法: {'全局' if args.algorithm == 'global' else '局部'} | 序列类型: {seq_type}")
        print("="*80)
        
        align_func = needleman_wunsch if args.algorithm == 'global' else smith_waterman
        results = batch_alignment(ref_seq, query_seqs, align_func, scoring)
        
        print(f"\n{'序号':<6}{'序列名称':<30}{'长度':<10}{'得分':<12}{'相似度(%)':<12}{'匹配数':<10}{'Gap数':<10}")
        print("-" * 100)
        
        for i, result in enumerate(results, 1):
            print(f"{i:<6}{result['name']:<30}{result['query_length']:<10}{result['score']:<12.1f}"
                  f"{result['similarity']:<12.2f}{result['matches']:<10}{result['gaps']:<10}")
        
        if args.output:
            export_batch_results(results, args.output, 
                                reference_name=ref_name,
                                algorithm=args.algorithm,
                                fmt=args.batch_format)
            print(f"\n结果已导出到: {args.output} (格式: {args.batch_format.upper()})")
        
        return
    
    seq1, seq2, name1, name2 = get_sequences(args)
    
    if not seq1 or not seq2:
        print("错误: 请提供两个序列进行比对", file=sys.stderr)
        print("使用 -h 查看帮助信息", file=sys.stderr)
        sys.exit(1)
    
    seq_type = args.seq_type
    if seq_type == 'auto':
        seq_type1 = detect_seq_type(seq1)
        seq_type2 = detect_seq_type(seq2)
        if seq_type1 != seq_type2:
            print(f"警告: 两个序列类型检测不一致 ({seq_type1} vs {seq_type2})，使用 {seq_type1}", file=sys.stderr)
        seq_type = seq_type1
    
    scoring = ScoringMatrix(
        seq_type=seq_type,
        match=args.match,
        mismatch=args.mismatch,
        gap=args.gap,
        matrix_name=args.matrix,
        matrix_file=args.matrix_file
    )
    
    print("\n" + "="*70)
    print("生物序列比对工具")
    print("="*70)
    print(f"序列1: {name1} (长度: {len(seq1)})")
    print(f"序列2: {name2} (长度: {len(seq2)})")
    print(f"算法: {'全局比对 (Needleman-Wunsch)' if args.algorithm == 'global' else '局部比对 (Smith-Waterman)'}")
    print(f"序列类型: {seq_type.upper()}")
    if scoring.matrix_name:
        print(f"打分矩阵: {scoring.matrix_name.upper()}")
    else:
        print(f"打分参数: match={args.match}, mismatch={args.mismatch}, gap={args.gap}")
    print("-"*70)
    
    align_func = needleman_wunsch if args.algorithm == 'global' else smith_waterman
    
    if args.reverse_complement:
        if seq_type != 'dna':
            print("警告: --reverse-complement 仅适用于DNA序列，将忽略此选项", file=sys.stderr)
            result = align_func(seq1, seq2, scoring)
            direction_info = ""
        else:
            print("\n正在比较正向和反向互补序列，选择最佳方向...", file=sys.stderr)
            rc_result = align_with_reverse_complement(seq1, seq2, scoring, algorithm=args.algorithm)
            result = rc_result['result']
            direction = rc_result['direction']
            forward_score = rc_result['forward_score']
            reverse_score = rc_result['reverse_score']
            
            if direction == 'reverse_complement':
                direction_info = f"\n方向: 使用序列2的反向互补 (正向得分: {forward_score:.1f}, 反向互补得分: {reverse_score:.1f})"
                name2 = f"{name2} (反向互补)"
            else:
                direction_info = f"\n方向: 使用原序列 (正向得分: {forward_score:.1f}, 反向互补得分: {reverse_score:.1f})"
    else:
        result = align_func(seq1, seq2, scoring)
        direction_info = ""
    
    print(f"\n比对得分: {result.score:.1f}")
    print(f"相似度: {result.similarity:.2f}%")
    print(f"匹配数: {result.matches} / {result.aligned_length}")
    print(f"错配数: {result.mismatches}")
    print(f"Gap数: {result.gaps}")
    if direction_info:
        print(direction_info)
    print("\n比对结果:")
    print("-"*70)
    
    print_alignment(result, name1=name1, name2=name2, line_width=args.line_width, 
                   use_color=not args.no_color)
    
    if args.output:
        export_alignment(result, args.output, args.format, name1=name1, name2=name2)
        print(f"\n结果已导出到: {args.output}")


if __name__ == '__main__':
    main()
