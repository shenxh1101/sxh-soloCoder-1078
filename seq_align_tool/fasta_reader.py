"""
FASTA文件读取模块
支持读取单个FASTA文件和目录中的所有FASTA文件
"""

import os
import re


def read_fasta(file_path):
    """
    读取FASTA格式的序列文件
    
    Args:
        file_path: FASTA文件路径
        
    Returns:
        列表，每个元素为 (name, sequence) 元组
    """
    sequences = []
    current_name = None
    current_seq = []
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_name is not None:
                        seq = ''.join(current_seq).upper().replace(' ', '').replace('\t', '')
                        sequences.append((current_name, seq))
                    current_name = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_name is not None:
                seq = ''.join(current_seq).upper().replace(' ', '').replace('\t', '')
                sequences.append((current_name, seq))
    except Exception as e:
        print(f"读取FASTA文件 {file_path} 时出错: {e}")
        return []
    
    return sequences


def read_fasta_directory(directory_path):
    """
    读取目录中所有FASTA格式的文件
    
    Args:
        directory_path: 目录路径
        
    Returns:
        列表，每个元素为 (name, sequence) 元组
    """
    sequences = []
    fasta_extensions = ['.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn']
    
    if not os.path.isdir(directory_path):
        print(f"错误: {directory_path} 不是有效的目录")
        return []
    
    for filename in sorted(os.listdir(directory_path)):
        file_path = os.path.join(directory_path, filename)
        if os.path.isfile(file_path):
            _, ext = os.path.splitext(filename)
            if ext.lower() in fasta_extensions:
                seqs = read_fasta(file_path)
                if seqs:
                    sequences.extend(seqs)
    
    return sequences


def validate_sequence(sequence, seq_type='dna'):
    """
    验证序列的有效性
    
    Args:
        sequence: 序列字符串
        seq_type: 序列类型 ('dna' 或 'protein')
        
    Returns:
        (is_valid, message) 元组
    """
    if not sequence:
        return False, "序列为空"
    
    if seq_type == 'dna':
        valid_chars = set('ATCGN')
        seq_upper = sequence.upper().replace('-', '')
        invalid_chars = set(seq_upper) - valid_chars
        if invalid_chars:
            return False, f"DNA序列包含无效字符: {', '.join(sorted(invalid_chars))}"
    elif seq_type == 'protein':
        valid_chars = set('ACDEFGHIKLMNPQRSTVWY*')
        seq_upper = sequence.upper().replace('-', '')
        invalid_chars = set(seq_upper) - valid_chars
        if invalid_chars:
            return False, f"蛋白质序列包含无效字符: {', '.join(sorted(invalid_chars))}"
    
    return True, "序列有效"


def write_fasta(sequences, output_path, line_width=80):
    """
    将序列写入FASTA文件
    
    Args:
        sequences: 列表，每个元素为 (name, sequence) 元组
        output_path: 输出文件路径
        line_width: 每行序列的最大宽度
    """
    with open(output_path, 'w', encoding='utf-8') as f:
        for name, seq in sequences:
            f.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')
