import argparse
import scrublet as scr
import numpy as np
import scipy.io
import matplotlib
import matplotlib.pyplot as plt
import os
import pandas as pd
import numba
import numba.typed

# 读取包含单细胞矩阵路径的文件列表
def read_matrix_paths(file_list):
    """
    读取包含单细胞矩阵路径的文本文件。

    参数：
        file_list (str): 包含矩阵路径的文本文件路径。

    返回：
        list: 包含矩阵路径的列表。
    """
    with open(file_list, 'r') as f:
        matrix_paths = [line.strip() for line in f.readlines()] # 读取每一行，去除首尾空白字符
    return matrix_paths

# 进行双细胞检测
def run_scrublet(matrix_paths, output_dir, doublet_score_threshold=None):
    """
    对每个单细胞矩阵进行双细胞检测。

    参数：
        matrix_paths (list): 包含矩阵路径的列表。
        output_dir (str): 输出结果的目录。
        doublet_score_threshold (float, 可选): 双细胞得分阈值。默认为None（自动）。
    """
    for matrix_path in matrix_paths:
        # 读取单细胞矩阵
        try:
            counts_matrix = scipy.io.mmread(matrix_path + '/matrix.mtx').T.tocsc() # 读取稀疏矩阵
        except FileNotFoundError:
            print(f"File not found: {matrix_path}/matrix.mtx")
            continue

        try:
            genes = np.array(scr.load_genes(matrix_path + '/features.tsv', delimiter='\t', column=1)) # 读取基因列表
        except FileNotFoundError:
            print(f"Genes file not found: {matrix_path}/genes.tsv")
            continue

        try:
            barcode = pd.read_csv(f"{matrix_path}/barcodes.tsv", sep="\t", header=None, names=["barcode"]) # 读取条形码
        except FileNotFoundError:
            print(f"Genes file not found: {matrix_path}/barcodes.tsv")
            continue
        

        # 提取样本名称
        parts = matrix_path.split('/')
        sample_name = parts[-1]  # 倒数第四个目录是样本名称
        
        print(f'Processing sample: {sample_name}')

        print(f'Counts matrix shape: {counts_matrix.shape[0]} rows, {counts_matrix.shape[1]} columns')
        print(f'Number of genes in gene list: {len(genes)}')
        print(f'Number of cell : {len(barcode)}')

        # 创建 Scrublet 实例
        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06, sim_doublet_ratio=2) # 初始化Scrublet对象

        # 运行双细胞检测
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                                    min_cells=3,
                                                                    min_gene_variability_pctl=85,
                                                                    n_prin_comps=30) # 运行双细胞检测

        # 使用 call_doublets 方法根据给定的双重细胞评分阈值进行预测
        if doublet_score_threshold is not None:
            predicted_doublets = scrub.call_doublets(threshold=doublet_score_threshold) # 使用指定阈值预测
        else:
            predicted_doublets = scrub.predicted_doublets_ # 使用自动阈值预测

        # 获取UMAP降维结果
        print('Running UMAP...')
        scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3)) # 运行UMAP降维
        print('Done.')

        # 绘制histogram
        scrub.plot_histogram() # 绘制双细胞得分直方图
        plt.savefig(f'{output_dir}/{sample_name}_doublet_histgram.png') # 保存直方图

        # 绘制双细胞预测的2D UMAP嵌入
        scrub.plot_embedding('UMAP', order_points=True) # 绘制UMAP图
        plt.savefig(f'{output_dir}/{sample_name}_doublet_predictions_umap.png') # 保存UMAP图

        # 保存双细胞检测结果
        results = pd.Series(predicted_doublets, name="scrublet_DropletType") # 创建双细胞预测结果Series
        scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores") # 创建双细胞得分Series
        dataframe = pd.concat([barcode, results, scores], axis=1) # 合并结果
        dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet") # 将True替换为doublet
        dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet") # 将False替换为singlet

        dataframe.to_csv(f'{output_dir}/{sample_name}_doublets.txt', sep="\t", index=False) # 保存结果到文件

        # Make summary of singlets and doublets and write to file
        summary = pd.DataFrame(dataframe.scrublet_DropletType.value_counts()) # 统计singlet和doublet数量
        summary.index.name = 'Classification'
        summary.reset_index(inplace=True)
        summary = summary.rename(columns={'scrublet_DropletType': 'DropletN'})
        summary.to_csv(f'{output_dir}/{sample_name}_doublet_summary.txt', sep="\t", index=False) # 保存统计结果

# 主函数
def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='Run Scrublet for doublet detection on single-cell matrices.')
    parser.add_argument('-s', '--single_cell_matrix_paths', required=True, help='Path to the text file containing a list of single-cell matrix paths')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory to save the output results')
    parser.add_argument('--doublet_score_threshold', type=float, default=None, help='Threshold for doublet score. Default is None (auto)')

    # 获取命令行参数
    args = parser.parse_args()

    matrix_list_file = args.single_cell_matrix_paths
    output_directory = args.output_dir
    doublet_score_threshold = args.doublet_score_threshold

    # 创建输出目录（如果不存在）
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # 读取矩阵路径并进行双细胞检测
    matrix_paths = read_matrix_paths(matrix_list_file)
    run_scrublet(matrix_paths, output_directory, doublet_score_threshold)

if __name__ == '__main__':
    main()
