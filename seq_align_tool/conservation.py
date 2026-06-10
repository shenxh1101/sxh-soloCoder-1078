"""
保守性分析模块
支持多序列保守性分析和ASCII保守度曲线图
"""

from typing import List, Dict, Tuple
from collections import Counter
from .alignment import needleman_wunsch
from .scoring import ScoringMatrix


def _align_multiple_sequences(sequences: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """
    渐进式多序列比对
    确保所有序列（包括master）对齐到同一组列，gap位置正确对齐
    """
    if len(sequences) < 2:
        return sequences
    
    seq_type = _detect_seq_type(sequences[0][1])
    scoring = ScoringMatrix(seq_type=seq_type, match=2, mismatch=-1, gap=-2)
    
    aligned = [(sequences[0][0], sequences[0][1])]
    
    for name, seq in sequences[1:]:
        current_master = aligned[0][1]
        result = needleman_wunsch(current_master, seq, scoring)
        
        aligned_master = result.seq1_aligned
        aligned_seq = result.seq2_aligned
        
        if len(aligned_master) != len(aligned[0][1]):
            aligned = _realign_all_sequences(aligned, aligned_master, result)
        
        aligned.append((name, aligned_seq))
    
    max_len = max(len(seq) for _, seq in aligned)
    aligned = [(name, seq.ljust(max_len, '-')) for name, seq in aligned]
    
    return aligned


def _realign_all_sequences(aligned, new_master, result):
    """
    当master序列因新gap插入而长度变化时，重新对齐所有已有序列
    """
    old_master = aligned[0][1]
    new_master = result.seq1_aligned
    
    mapping = []
    old_idx = 0
    for c in new_master:
        if c == '-':
            mapping.append(None)
        else:
            mapping.append(old_idx)
            old_idx += 1
    
    new_aligned = []
    for name, seq in aligned:
        new_seq = []
        seq_idx = 0
        for pos in mapping:
            if pos is None:
                new_seq.append('-')
            else:
                if seq_idx < len(seq):
                    new_seq.append(seq[seq_idx])
                    seq_idx += 1
                else:
                    new_seq.append('-')
        new_aligned.append((name, ''.join(new_seq)))
    
    new_aligned[0] = (aligned[0][0], new_master)
    
    return new_aligned


def _detect_seq_type(seq: str) -> str:
    """检测序列类型"""
    dna_chars = set('ATCGN')
    seq_upper = seq.upper().replace('-', '')
    if all(c in dna_chars for c in seq_upper):
        return 'dna'
    return 'protein'


def calculate_conservation(aligned_sequences: List[Tuple[str, str]]) -> Dict:
    """
    计算多序列比对每个位置的保守度
    所有位置（包括gap列）都计入统计
    
    Args:
        aligned_sequences: 比对后的序列列表，每个元素为 (name, aligned_seq) 元组
        
    Returns:
        包含保守性信息的字典
    """
    if len(aligned_sequences) < 2:
        raise ValueError("保守性分析需要至少2个序列")
    
    max_len = max(len(seq) for _, seq in aligned_sequences)
    num_sequences = len(aligned_sequences)
    
    aligned_padded = []
    for name, seq in aligned_sequences:
        padded = seq.ljust(max_len, '-')
        aligned_padded.append((name, padded))
    
    conservation_scores = []
    conservation_levels = []
    gap_counts = []
    
    for pos in range(max_len):
        residues = []
        for _, seq in aligned_padded:
            residues.append(seq[pos])
        
        count = Counter(residues)
        gap_count = count.get('-', 0)
        gap_counts.append(gap_count)
        
        non_gap_residues = [r for r in residues if r != '-']
        
        if gap_count == num_sequences:
            conservation_scores.append(0.0)
            conservation_levels.append('all_gap')
            continue
        
        if not non_gap_residues:
            conservation_scores.append(0.0)
            conservation_levels.append('all_gap')
            continue
        
        non_gap_count = Counter(non_gap_residues)
        most_common_residue, most_common_count = non_gap_count.most_common(1)[0]
        
        conservation = (most_common_count / len(non_gap_residues)) * 100
        conservation_scores.append(conservation)
        
        if conservation == 100:
            conservation_levels.append('conserved')
        elif conservation >= 80:
            conservation_levels.append('high')
        elif conservation >= 50:
            conservation_levels.append('moderate')
        else:
            conservation_levels.append('low')
    
    valid_scores = [s for s, l in zip(conservation_scores, conservation_levels) if l != 'all_gap']
    avg_conservation = sum(valid_scores) / len(valid_scores) if valid_scores else 0
    conserved_sites = sum(1 for s in conservation_scores if s == 100)
    all_gap_sites = sum(1 for l in conservation_levels if l == 'all_gap')
    
    return {
        'positions': list(range(1, max_len + 1)),
        'conservation_scores': conservation_scores,
        'conservation_levels': conservation_levels,
        'gap_counts': gap_counts,
        'avg_conservation': avg_conservation,
        'conserved_sites': conserved_sites,
        'all_gap_sites': all_gap_sites,
        'total_sites': max_len,
        'num_sequences': num_sequences,
        'aligned_sequences': aligned_padded
    }


def generate_ascii_plot(conservation_scores: List[float], plot_width: int = 80, 
                       height: int = 10) -> str:
    """
    生成保守度的ASCII曲线图
    
    Args:
        conservation_scores: 保守度分数列表
        plot_width: 图表宽度（字符数）
        height: 图表高度（行数）
        
    Returns:
        ASCII图表字符串
    """
    if not conservation_scores:
        return "无数据"
    
    n = len(conservation_scores)
    
    if n > plot_width:
        bin_size = n // plot_width
        binned_scores = []
        for i in range(0, n, bin_size):
            chunk = conservation_scores[i:i+bin_size]
            if chunk:
                binned_scores.append(sum(chunk) / len(chunk))
        conservation_scores = binned_scores
        n = len(conservation_scores)
    
    plot_lines = []
    
    for h in range(height, -1, -1):
        threshold = h * (100.0 / height)
        line = '|'
        
        for score in conservation_scores:
            if score >= threshold:
                if score >= 100:
                    line += '*'
                elif score >= 80:
                    line += '#'
                elif score >= 50:
                    line += '+'
                else:
                    line += '.'
            else:
                line += ' '
        
        if h % 2 == 0:
            line += f' {int(threshold):3d}%'
        plot_lines.append(line)
    
    x_axis = '+' + '-' * len(conservation_scores) + '+'
    plot_lines.append(x_axis)
    
    if n > 0:
        label_line = '|'
        step = max(1, n // 5)
        for i in range(n):
            if i % step == 0:
                label_line += str((i + 1) % 10)
            else:
                label_line += ' '
        label_line += '| Position'
        plot_lines.append(label_line)
    
    legend = []
    legend.append("")
    legend.append("图例: * = 100% (完全保守), # = >=80%, + = >=50%, . = <50%")
    legend.append(f"平均保守度: {sum(conservation_scores)/len(conservation_scores):.1f}%")
    legend.append(f"完全保守位点: {sum(1 for s in conservation_scores if s >= 100)} / {len(conservation_scores)}")
    
    plot_lines.extend(legend)
    
    return '\n'.join(plot_lines)


def conservation_analysis(sequences: List[Tuple[str, str]], plot_width: int = 80) -> Dict:
    """
    执行完整的保守性分析
    
    Args:
        sequences: 序列列表，每个元素为 (name, sequence) 元组
        plot_width: ASCII曲线图宽度
        
    Returns:
        包含保守性分析结果的字典
    """
    if len(sequences) < 2:
        raise ValueError("保守性分析需要至少2个序列")
    
    aligned = _align_multiple_sequences(sequences)
    
    result = calculate_conservation(aligned)
    
    ascii_plot = generate_ascii_plot(result['conservation_scores'], plot_width=plot_width)
    result['ascii_plot'] = ascii_plot
    result['aligned_sequences'] = aligned
    
    return result


def print_conservation_result(result: Dict) -> None:
    """
    打印保守性分析结果
    """
    print("\n" + "="*70)
    print("多序列保守性分析结果")
    print("="*70)
    print(f"\n序列数量: {result['num_sequences']}")
    print(f"比对总长度: {result['total_sites']} 个位置")
    print(f"完全保守位点: {result['conserved_sites']} / {result['total_sites']} "
          f"({result['conserved_sites']/result['total_sites']*100:.1f}%)")
    print(f"全gap位点: {result['all_gap_sites']} / {result['total_sites']}")
    print(f"有效保守位点（排除全gap）: {result['total_sites'] - result['all_gap_sites']}")
    print(f"平均保守度（仅有效位点）: {result['avg_conservation']:.2f}%")
    
    if 'aligned_sequences' in result:
        print(f"\n多序列比对结果（前60个碱基）:")
        print("-" * 70)
        for name, seq in result['aligned_sequences']:
            display_seq = seq[:60] + ('...' if len(seq) > 60 else '')
            print(f"{name[:25]:<25}: {display_seq}")
    
    print(f"\n保守度曲线图:")
    print("-" * 70)
    print(result['ascii_plot'])
    print("=" * 70)
