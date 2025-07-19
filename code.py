import requests
from Bio import Entrez
import pandas as pd
import xml.etree.ElementTree as ET
from typing import List, Dict
import time
import urllib.request
import urllib.error
import http.client
import ssl

class FluDataCrawler:
    def __init__(self, email: str):
        """
        初始化爬虫，使用您的电子邮件（NCBI E-utilities 要求）
        """
        Entrez.email = email
        Entrez.api_key = None  # 可选：如果有 NCBI API 密钥，请在此设置
        self.base_url_ird = "https://www.fludb.org/brc/api/"
        self.headers = {'User-Agent': 'FluDataCrawler/1.0 (Python; +your.email@example.com)'}

    def search_ncbi_flu(self, query: str, retmax: int = 500) -> List[str]:
        """
        搜索 NCBI 流感病毒数据库
        Args:
            query: 搜索查询字符串
            retmax: 最大返回结果数（默认：500）
        """
        print(f"正在搜索最多 {retmax} 个序列...")
        for attempt in range(3):  # 重试 3 次
            try:
                handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax, usehistory="y")
                record = Entrez.read(handle)
                handle.close()
                return record["IdList"]
            except (urllib.error.URLError, http.client.HTTPException) as e:
                print(f"搜索失败，重试 {attempt + 1}/3: {str(e)}")
                time.sleep(2 ** attempt)  # 指数退避
            except Exception as e:
                print(f"搜索时发生未知错误: {str(e)}")
                return []
        print("搜索失败，超过最大重试次数。")
        return []

    def fetch_sequence_data(self, id_list: List[str]) -> List[Dict]:
        """
        获取给定 ID 的序列和注释数据
        """
        sequences = []
        batch_size = 50  # 减小批处理大小，减少服务器压力

        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i + batch_size]
            for attempt in range(3):  # 重试 3 次
                try:
                    handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="gb", retmode="xml")
                    records = Entrez.parse(handle)

                    for record in records:
                        seq_data = {
                            'accession': record.get('GBSeq_accession-version', ''),
                            'sequence': record.get('GBSeq_sequence', ''),
                            'strain': '',
                            'host': '',
                            'collection_date': '',
                            'country': '',
                            'virulence_info': ''
                        }

                        # 提取特征和限定词
                        for feature in record.get('GBSeq_feature-table', []):
                            if feature.get('GBFeature_key') == 'source':
                                for qualifier in feature.get('GBFeature_quals', []):
                                    if qualifier.get('GBQualifier_name') == 'strain':
                                        seq_data['strain'] = qualifier.get('GBQualifier_value', '')
                                    elif qualifier.get('GBQualifier_name') == 'host':
                                        seq_data['host'] = qualifier.get('GBQualifier_value', '')
                                    elif qualifier.get('GBQualifier_name') == 'country':
                                        seq_data['country'] = qualifier.get('GBQualifier_value', '')
                                    elif qualifier.get('GBQualifier_name') == 'collection_date':
                                        seq_data['collection_date'] = qualifier.get('GBQualifier_value', '')

                        # 在注释和备注中查找毒力信息
                        for feature in record.get('GBSeq_feature-table', []):
                            if 'GBFeature_quals' in feature:
                                for qualifier in feature['GBFeature_quals']:
                                    if qualifier.get('GBQualifier_name') in ['note', 'comment']:
                                        value = qualifier.get('GBQualifier_value', '').lower()
                                        if any(term in value for term in ['virulen', 'pathogenic', 'pathogenicity']):
                                            seq_data['virulence_info'] = qualifier.get('GBQualifier_value', '')

                        sequences.append(seq_data)

                    handle.close()
                    print(f"成功处理批次 {i // batch_size + 1}")
                    time.sleep(1)  # 对 NCBI 服务器友好
                    break
                except Exception as e:
                    print(f"处理批次 {i // batch_size + 1} 失败，重试 {attempt + 1}/3: {str(e)}")
                    time.sleep(2 ** attempt)
                    continue
            else:
                print(f"批次 {i // batch_size + 1} 失败，跳过。")
                continue

        return sequences

    def search_ird(self, params: Dict) -> List[Dict]:
        """
        搜索流感研究数据库（IRD）
        注意：这是一个基本实现，可能需要根据 IRD 的 API 进行调整
        """
        try:
            response = requests.get(f"{self.base_url_ird}/search", params=params, headers=self.headers, verify=False)
            if response.status_code == 200:
                return response.json()
            else:
                print(f"访问 IRD 失败: {response.status_code}")
                return []
        except Exception as e:
            print(f"访问 IRD 时发生错误: {str(e)}")
            return []

def main():
    # 初始化爬虫，使用您的电子邮件
    crawler = FluDataCrawler("your.email@example.com")

    # 搜索流感 A 病毒序列，包含毒力信息
    query = "(Influenza A virus[Organism]) AND (virulence[All Fields] OR pathogenicity[All Fields])"
    id_list = crawler.search_ncbi_flu(query)

    if id_list:
        print(f"找到 {len(id_list)} 个序列，开始下载...")
        # 获取每个序列的详细信息
        sequences = crawler.fetch_sequence_data(id_list)

        if not sequences:
            print("未成功下载任何序列，请检查上述错误。")
            return

        # 转换为 DataFrame 以便于处理
        df = pd.DataFrame(sequences)

        # 保存完整信息到 CSV
        df.to_csv("flu_sequences_with_virulence.csv", index=False)
        print(f"成功保存 {len(sequences)} 个序列到 flu_sequences_with_virulence.csv")

        # 以 FASTA 格式保存序列
        with open("flu_sequences.fasta", "w") as f:
            for seq in sequences:
                f.write(f">{seq['accession']}|{seq['strain']}|{seq[' host']}|{seq['collection_date']}\n")
                f.write(f"{seq['sequence']}\n")
        print("成功保存序列到 flu_sequences.fasta")

        # 创建并保存 ID-毒力注释 CSV
        virulence_df = pd.DataFrame({
            'sequence_id': df['accession'],
            'virulence_annotation': df['virulence_info']
        })
        # 仅保留包含毒力信息的条目
        virulence_df = virulence_df[virulence_df['virulence_annotation'].str.len() > 0]
        virulence_df.to_csv("virulence_annotations.csv", index=False)
        print(f"成功保存 {len(virulence_df)} 个毒力注释到 virulence_annotations.csv")
    else:
        print("未找到任何序列，请检查您的搜索查询。")

if __name__ == "__main__":
    main()
