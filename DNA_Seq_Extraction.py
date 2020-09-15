import time

from Bio.Seq import Seq
import re
import requests


# 3.15
def getHTMLText(url):
    """
    请求接面
    :param url: 拼接的url，利用接口直接返回数据
    :return: html
    """
    try:
        header = {
            'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/78.0.3904.108 Safari/537.36'
        }
        response = requests.get(url, headers=header, timeout=30)
        # 异常判断
        response.raise_for_status()
        # 获取编码格式
        response.encoding = response.apparent_encoding
        return response.text
    except:
        return ""


def fillDNAlist(html):
    """
    解析html页面，获取需要bp
    :param html:
    :return:
    """
    try:
        res = re.findall(r'<DNA\s+length="\d+">\s+(.*\s+\w+)\s+</DNA>', html)[0]
        res = res.replace("\n", "").upper()
    except Exception as e:
        res = ""
    return res


def variant_change(dna_seq, chrom, pos, normal, variant):
    """
    :param dna_seq: 正确的基因序列
    :param chrom: 染色体号
    :param pos: 染色体位置
    :param normal: 正常基因
    :param variant: 突变基因
    :return:
    """
    if dna_seq == "":
        pass
    else:
        print("var in res DNA Seq", dna_seq)
        # 突变基因两端基因片段
        left = dna_seq[:25]
        right = dna_seq[26:]

        # 利用Seq对象处理基因序列，基因互补配对【逆序】
        my_seq = Seq(right)
        left_3 = my_seq.reverse_complement()
        my_seq = Seq(left)
        right_3 = my_seq.reverse_complement()

        with open(f"Y088N_data.fasta", "a") as f_w:
            print("开始写入---------")
            row_one = f">{chrom}_{pos}_{normal}_{variant}|+|ref" + "\n"
            row_two = dna_seq + "\n"
            row_there = f">{chrom}_{pos}_{normal}_{variant}|+|alt" + "\n"
            row_four = left + variant + right + "\n"
            row_five = f">{chrom}_{pos}_{normal}_{variant}|-|ref" + "\n"
            row_six = left_3 + Seq(normal).complement() + right_3 + "\n"
            row_seven = f">{chrom}_{pos}_{normal}_{variant}|-|alt" + "\n"
            row_eight = left_3 + Seq(variant).complement() + right_3 + "\n"

            res_data = [str(row_one), str(row_two), str(row_there), str(row_four), str(row_five), str(row_six),
                        str(row_seven), str(row_eight)]
            f_w.writelines(res_data)
            f_w.writelines("\n")
            f_w.close()
        print("写入结束--------")


if __name__ == '__main__':
    start_time = time.time()
    hg19 = "hg19"
    with open("Y088N_data.vcf") as file:
        data = file.readlines()
        for item in data:
            row = item.split()
            # 注意row[0]的格式，具体问题具体分析
            # 如果染色体仅是单独的数字
            chrom = row[0]
            # 如果染色体是chr+数字
            # chrom = row[0][3:]
            pos = int(row[1])
            start = pos - 25
            end = pos + 25
            normal = row[3]
            variant = row[4]

            # 请求连接拼接
            # http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr17:7675591,7676591
            url = f"http://genome.ucsc.edu/cgi-bin/das/{hg19}/dna?segment=chr{chrom}:{start},{end}"
            print(url)

            html = getHTMLText(url)
            dna_seq = fillDNAlist(html)
            variant_change(dna_seq, chrom, pos, normal, variant)

    end_time = time.time()

    print("run final！！！")
    print("cost time:", end_time - start_time)
