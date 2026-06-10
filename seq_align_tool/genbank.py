"""
GenBank序列获取模块
支持通过GenBank accession号自动获取序列（需联网）
"""

import sys
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET


def fetch_genbank(accession: str) -> tuple:
    """
    从NCBI GenBank获取序列
    
    Args:
        accession: GenBank accession号 (如: NM_001301717, U49845)
        
    Returns:
        (name, sequence) 元组，如果获取失败返回 (None, None)
    """
    try:
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        
        params = {
            'db': 'nucleotide',
            'id': accession,
            'rettype': 'fasta',
            'retmode': 'text',
            'tool': 'seq_align_tool',
            'email': 'user@example.com'
        }
        
        url = f"{base_url}?{urllib.parse.urlencode(params)}"
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (compatible; SeqAlignTool/1.0)'
        }
        
        req = urllib.request.Request(url, headers=headers)
        
        with urllib.request.urlopen(req, timeout=30) as response:
            if response.status != 200:
                print(f"错误: HTTP状态码 {response.status}", file=sys.stderr)
                return None, None
            
            data = response.read().decode('utf-8')
            
            if not data or data.startswith('Error') or 'Error' in data[:100]:
                print(f"错误: 无法获取序列 {accession}", file=sys.stderr)
                return None, None
            
            lines = data.strip().split('\n')
            if not lines or not lines[0].startswith('>'):
                print(f"错误: 返回的数据不是有效的FASTA格式", file=sys.stderr)
                return None, None
            
            name = lines[0][1:].strip()
            sequence = ''.join(lines[1:]).upper().replace(' ', '').replace('\t', '')
            
            if not sequence:
                print(f"错误: 序列为空", file=sys.stderr)
                return None, None
            
            return name, sequence
            
    except urllib.error.HTTPError as e:
        print(f"HTTP错误: {e.code} - {e.reason}", file=sys.stderr)
        return None, None
    except urllib.error.URLError as e:
        print(f"URL错误: {e.reason}", file=sys.stderr)
        print("请检查网络连接或稍后重试", file=sys.stderr)
        return None, None
    except Exception as e:
        print(f"获取GenBank序列时发生错误: {e}", file=sys.stderr)
        return None, None


def fetch_genbank_protein(accession: str) -> tuple:
    """
    从NCBI GenBank获取蛋白质序列
    
    Args:
        accession: GenBank蛋白质accession号
        
    Returns:
        (name, sequence) 元组，如果获取失败返回 (None, None)
    """
    try:
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        
        params = {
            'db': 'protein',
            'id': accession,
            'rettype': 'fasta',
            'retmode': 'text',
            'tool': 'seq_align_tool',
            'email': 'user@example.com'
        }
        
        url = f"{base_url}?{urllib.parse.urlencode(params)}"
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (compatible; SeqAlignTool/1.0)'
        }
        
        req = urllib.request.Request(url, headers=headers)
        
        with urllib.request.urlopen(req, timeout=30) as response:
            if response.status != 200:
                print(f"错误: HTTP状态码 {response.status}", file=sys.stderr)
                return None, None
            
            data = response.read().decode('utf-8')
            
            if not data or data.startswith('Error') or 'Error' in data[:100]:
                print(f"错误: 无法获取序列 {accession}", file=sys.stderr)
                return None, None
            
            lines = data.strip().split('\n')
            if not lines or not lines[0].startswith('>'):
                print(f"错误: 返回的数据不是有效的FASTA格式", file=sys.stderr)
                return None, None
            
            name = lines[0][1:].strip()
            sequence = ''.join(lines[1:]).upper().replace(' ', '').replace('\t', '')
            
            if not sequence:
                print(f"错误: 序列为空", file=sys.stderr)
                return None, None
            
            return name, sequence
            
    except Exception as e:
        print(f"获取GenBank蛋白质序列时发生错误: {e}", file=sys.stderr)
        return None, None


def search_genbank(term: str, db: str = 'nucleotide', max_results: int = 10) -> list:
    """
    在GenBank中搜索序列
    
    Args:
        term: 搜索关键词
        db: 数据库 ('nucleotide' 或 'protein')
        max_results: 最大返回结果数
        
    Returns:
        包含accession号和描述的列表
    """
    try:
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        
        params = {
            'db': db,
            'term': term,
            'retmax': max_results,
            'retmode': 'xml',
            'tool': 'seq_align_tool',
            'email': 'user@example.com'
        }
        
        url = f"{search_url}?{urllib.parse.urlencode(params)}"
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (compatible; SeqAlignTool/1.0)'
        }
        
        req = urllib.request.Request(url, headers=headers)
        
        with urllib.request.urlopen(req, timeout=30) as response:
            if response.status != 200:
                return []
            
            data = response.read().decode('utf-8')
            root = ET.fromstring(data)
            
            ids = [id_elem.text for id_elem in root.findall('.//Id')]
            
            if not ids:
                return []
            
            fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            params = {
                'db': db,
                'id': ','.join(ids),
                'rettype': 'fasta',
                'retmode': 'text',
                'tool': 'seq_align_tool',
                'email': 'user@example.com'
            }
            
            url = f"{fetch_url}?{urllib.parse.urlencode(params)}"
            req = urllib.request.Request(url, headers=headers)
            
            with urllib.request.urlopen(req, timeout=30) as response2:
                data2 = response2.read().decode('utf-8')
                
                results = []
                entries = data2.split('>')[1:]
                
                for entry in entries[:max_results]:
                    lines = entry.split('\n')
                    if lines:
                        header = lines[0].strip()
                        accession = header.split()[0] if header else ''
                        description = header[len(accession):].strip() if accession else header
                        results.append({
                            'accession': accession,
                            'description': description
                        })
                
                return results
            
    except Exception as e:
        print(f"搜索GenBank时发生错误: {e}", file=sys.stderr)
        return []
