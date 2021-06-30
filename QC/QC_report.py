import os, sys, argparse
from fpdf import FPDF
import pandas as pd

def parse_args() -> argparse.Namespace:
    """
    Returns:
        arguments
    """
    parser = argparse.ArgumentParser()

    files = parser.add_argument_group('File Arguments')
    files.add_argument('--first_log_file', required=True)
    files.add_argument('--last_log_file', required=True)
    files.add_argument('--fail_sex_file', required=True)
    files.add_argument('--fail_het_file', required=True)
    files.add_argument('--fail_miss_file', required=True)
    files.add_argument('--fail_ibd_file', required=True)
    files.add_argument('--fail_ps_file', required=True)

    figs = parser.add_argument_group('Figure Arguments')
    figs.add_argument('--fig_het_vs_imiss', required=True)
    figs.add_argument('--fig_hwe', required=True)
    figs.add_argument('--fig_maf', required=True)
    figs.add_argument('--fig_lmiss', required=True)
    figs.add_argument('--fig_population', required=True)

    params = parser.add_argument_group('Parameters')
    params.add_argument("--heter", type=float, required=True)
    params.add_argument("--ind_miss", type=float, required=True)
    params.add_argument("--geno_miss", type=float, required=True)
    params.add_argument("--ibd", type=float, required=True)
    params.add_argument("--maf", type=float, required=True)
    params.add_argument("--hwe", type=float, required=True)
    params.add_argument("--prune_window", type=float, required=True)
    params.add_argument("--prune_step", type=float, required=True)
    params.add_argument("--prune_threshold", type=float, required=True)
    params.add_argument("--population", type=str, required=False, default="")

    parser.add_argument("--logo_file", type=str, required=True)
    parser.add_argument("--work_dir", type=str, required=True)
    parser.add_argument("--basename", type=str, required=True)

    args = parser.parse_args()

    return args


class PDF(FPDF):
    def page_header(self, Title):
        # Arial bold 15
        self.set_font('Times', 'B', 18)
        # Calculate width of title and position
        w = self.get_string_width(Title) + 6
        self.set_xy((210 - w) / 2, 15)
        # Colors of frame, background and text
        self.set_draw_color(256, 256, 256)
        self.set_fill_color(256, 256, 256)
        self.set_text_color(22, 143, 153)
        # Thickness of frame (1 mm)
        self.set_line_width(1)
        # Title
        self.cell(w, 9, Title, 1, 1, 'C', 1)
        # Line break
        self.ln(10)
        
    def logo(self, fig):
        self.image(fig, x = 5, y = 5, w = 50, h = 12)
        
    def footer(self):
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        # Arial italic 8
        self.set_font('Times', 'I', 8)
        # Text color in gray
        self.set_text_color(51, 51, 51)
        # Page number
        self.cell(0, 10, 'Page ' + str(self.page_no()), 0, 0, 'C')

    def chapter_title(self, num, label):
        # Arial 12
        self.set_font('Times', 'B', 13)
        # Background color
        self.set_fill_color(22, 143, 153)
        self.set_text_color(256, 256, 256)
        # Title
        self.cell(0, 6, ' Section %s : %s' % (num, label), 0, 1, 'L', 1)
        # Line break
        self.ln(4)

    def chapter_body(self, txt):
        # Times 12
        self.set_font('Times', '', 12)
        # Output justified text
        self.multi_cell(0, 5, txt)

    def content(self, text_format_list):
        for text, form in text_format_list:
            if form == 'bold':
                self.set_font('Times', 'B', 12)
                self.write(5, text)
            else:
                self.set_font('Times', '', 12)
                self.write(5, text)


def main():
    # Arguments
    args = parse_args()
    heter = args.heter
    ind_miss = args.ind_miss
    geno_miss = args.geno_miss
    ibd = args.ibd
    maf = args.maf
    hwe = args.hwe
    prune_window = args.prune_window
    prune_threshold = args.prune_threshold
    population = args.population
    basename = args.basename
    logo_file = args.logo_file

    os.chdir(args.work_dir)
    
    pdf = PDF()
    Title = "GWAS QC Report"
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(logo_file)
    pdf.set_line_width(0.5)
    
    ############## Introduction ################################
    # title
    pdf.chapter_title("1", "Introduction")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # get the number of samples and SNPs
    with open(args.first_log_file, 'r') as f:
        for line in f.readlines():
            if 'males' in line:
                end_point = line.index('loaded') - 1
                sample_info = line[:end_point]
                sample_info = sample_info.replace('people', 'samples')
                sample_num = int(line.split()[0])
            elif 'variants loaded' in line:
                nums = [int(s) for s in line.split() if s.isdigit()]
                snp_num = nums[0]
            elif 'variants removed due to missing genotype data' in line:
                nums = [int(s) for s in line.split() if s.isdigit()]
                snp_fail_miss_num = nums[0]
            elif 'variants removed due to Hardy-Weinberg exact test' in line:
                nums = [int(s) for s in line.split() if s.isdigit()]
                snp_fail_hwe_num = nums[0]
            elif 'variants removed due to minor allele threshold(s)' in line:
                nums = [int(s) for s in line.split() if s.isdigit()]
                snp_fail_maf_num = nums[0]

    # text and format
    text_list = [
        ('This report encompasses the quality control (QC) summary for the ', 'normal'),
        (basename, 'bold'),
        (' data. A total of ', 'normal'),
        (sample_info, 'bold'),
        (' were genotyped for ', 'normal'),
        (str(snp_num), 'bold'),
        (' SNPs. Quality control was performed across samples, across SNPs and on physical locations. The results and the filtering criteria are discussed below.', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)

    ################## cutoffs setting ########################### 
    # title
    pdf.chapter_title("2", "User-defined/Default cutoffs")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)
    
    # text and format
    text_list = [
        ('-- Individual Missingness Cutoff:', 'normal'),
        ('{0} (upper bound)'.format(ind_miss), 'bold'),
        ('-- Heterozygosity Cutoff:', 'normal'),
        ('{0} SD (standard deviation)'.format(heter), 'bold'),
        ('-- Identify By Decent Cutoff:', 'normal'),
        ('{0} (upper bound)'.format(ibd), 'bold'),
        ('-- Minor Allele Frequency Cutoff:', 'normal'),
        ('{0} (lower bound)'.format(maf), 'bold'),
        ('-- Genotype Missingness Cutoff:', 'normal'),
        ('{0} (upper bound)'.format(geno_miss), 'bold'),
        ('-- Hardy Weinberg Equilibrium Cutoff:', 'normal'),
        ('{0} (lower bound)'.format(hwe), 'bold'),
    ]

    # content
    pdf.chapter_body("QC was processed using the following cutoffs:")
    pdf.ln(2)
    for text, form in text_list:
        if form == 'bold':
            pdf.set_font('Times', 'B', 12)
            pdf.cell(0, 5, text)
            pdf.ln(5)
        else:
            pdf.cell(7)
            pdf.set_font('Times', '', 12)
            pdf.cell(70, 5, text)
    pdf.ln(10)


    ################## Sample summary #################################
    # title
    pdf.chapter_title("3.1", "Summary of Individual level QC")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # get the number of failed samples
    with open(args.last_log_file, 'r') as f:
        for line in f.readlines():
            if 'Among remaining phenotypes' in line:
                nums = [int(s) for s in line.split() if s.isdigit()]
                case_num, control_num = nums[0], nums[1]
            elif 'variants and ' in line:
                nums = [int(s) for s in line.split() if s.isdigit()]
                snp_left_num = nums[0]
                control_num = nums[1]
                case_num = 0

    fail_sample_num = sample_num - case_num - control_num

    # text and format
    text_list = [
        ('{0}'.format(case_num + control_num), 'bold'),
        (' people pass filters and QC ({0} cases and {1} controls. '.format(case_num, control_num), 'normal'),
        ('{0}'.format(fail_sample_num), 'bold'),
        (' samples failed the following QC procedures and were excluded from the GWAS analysis. ', 'normal'),
        ('The table below shows these samples and the failure procedure for each sample.', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)

    # get the FID and IID of the failed samples
    fail_file_dict = {
        'Gender': args.fail_sex_file,
        'Missingness': args.fail_miss_file,
        'Heterozygosity': args.fail_het_file,
        'IBD': args.fail_ibd_file,
        'Population': args.fail_ps_file
    }
    fail_num_dict = {
        'Gender': 0,
        'Missingness': 0,
        'Heterozygosity': 0,
        'IBD': 0,
        'Population': 0
    }
    
    fail_df = pd.DataFrame(columns=['FID', 'IID'])
    for col, f in fail_file_dict.items():
        if os.path.getsize(f) == 0:
            fail_df[col] = 0
        else:
            temp_df = pd.read_csv(f, sep='\s+', header=None)
            temp_df = temp_df.rename(columns={0: 'FID', 1:'IID'})
            temp_df = temp_df[['FID', 'IID']]
            temp_df[col] = 1
            fail_df = fail_df.merge(temp_df, on=['FID', 'IID'], how='outer')
            fail_num_dict[col] = temp_df.shape[0]
    fail_df.fillna(0, inplace=True)
    

    # table
    if fail_sample_num != 0:
        ## column width
        inter_space = 3
        table_width = 200
        text_size = 13
        while table_width > 190:
            text_size -= 1
            pdf.set_font('Times', '', text_size)
            col_width_dict = {col: pdf.get_string_width(col) for col in fail_df.columns}
            for idx, row in fail_df.iterrows():
                for col in ['FID', 'IID']:
                    if pdf.get_string_width(str(row[col])) > col_width_dict[col]:
                        col_width_dict[col] = pdf.get_string_width(str(row[col]))
            table_width = sum(col_width_dict.values()) + inter_space * len(fail_df.columns) * 2
        
        ## header
        pdf.set_font('Times', 'B', text_size)
        pdf.cell((190 - table_width)/2)
        pdf.set_text_color(22, 143, 153)
        pdf.set_draw_color(22, 143, 153)
        for col in fail_df.columns:
            if (col == 'FID') or (col == 'IID'):
                pdf.cell(col_width_dict[col] + inter_space*2, 5, txt=col, align='C', border="B")
            else:
                pdf.cell(col_width_dict[col] + inter_space*2, 5, txt=col, align='C', border="LB")
        pdf.ln(6)

        ## row
        pdf.set_font('Times', '', text_size)
        pdf.set_text_color(51, 51, 51)
        for idx, row in fail_df.iterrows():
            pdf.cell((190 - table_width)/2)
            for col in fail_df.columns:
                if (col == 'FID') or (col == 'IID'):
                    pdf.set_text_color(51, 51, 51)
                    pdf.cell(col_width_dict[col] + inter_space*2, 5, txt=str(row[col]), align='C')
                else:
                    if row[col] == 1:
                        pdf.set_text_color(255, 0, 0)
                        pdf.cell(col_width_dict[col] + inter_space*2, 5, txt='X', align='C', border='L')
                    else:
                        pdf.cell(col_width_dict[col] + inter_space*2, 5, txt='', align='C', border='L')
            pdf.ln(6)
    pdf.ln(6)


    ################## SNPs summary ################################
    # title
    pdf.chapter_title("3.2", "Summary of SNPs level QC") 
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # text and format
    text_list = [
        (str(snp_num - snp_left_num), 'bold'),
        (' out of ', 'normal'),
        (str(snp_num), 'bold'),
        (' were removed from SNPs QC procedures including Minor Allele Frequency (MAF), Missingness of SNPs, Hardy Weinberg Equilibrium (HWE), and abnormal chromosome (including sex chromosome). ', 'normal'),
        ('The table below shows the number of excluded SNPs in each step.', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)
    
    # table
    ## column width
    header = ['', 'MAF', 'Missingness', 'HWE', 'Abnormal Chr', 'Total']
    snp_fail_chr_num = snp_num - snp_fail_maf_num - snp_fail_miss_num - snp_fail_hwe_num - snp_left_num
    row = ['number of SNPs', str(snp_fail_maf_num), str(snp_fail_miss_num), str(snp_fail_hwe_num), str(snp_fail_chr_num), str(snp_num)]
    col_width = [max(pdf.get_string_width(str(header[i])), pdf.get_string_width(str(row[i]))) for i in range(len(header))]
    
    ## header
    table_width = sum(col_width) + inter_space * len(header) * 2
    text_size = 12
    pdf.set_font('Times', 'B', text_size)
    pdf.cell((190 - table_width)/2)
    pdf.set_draw_color(22, 143, 153)
    pdf.set_text_color(22, 143, 153)
    for i in range(len(header)):
        if i == 0:
            pdf.cell(col_width[i] + 2*inter_space, 5, txt=header[i], align='C', border='B')
        else:
            pdf.cell(col_width[i] + 2*inter_space, 5, txt=header[i], align='C', border='LB')
    pdf.ln(6)

    ## row
    pdf.cell((190 - table_width)/2)
    pdf.set_text_color(51, 51, 51)
    for i in range(len(row)):
        if i == 0:
            pdf.set_font('Times', 'B', text_size)
            pdf.cell(col_width[i] + 2*inter_space, 5, txt=row[i], align='C', border='')
        else:
            pdf.set_font('Times', '', text_size)
            pdf.cell(col_width[i] + 2*inter_space, 5, txt=row[i], align='C', border='L')
    pdf.ln(6)

    ################## Sample QC #################################
    Title = "Individual level QC"
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(logo_file)

    # Gender
    ## title
    pdf.chapter_title("4.1", "Individual level QC --- Gender Discrepancy")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    ## text and format
    text_list = [
        ('Checks for gender discrepancies between the individuals recorded in the dataset and their gender based on X chromosome heterozygosity/homozygosity rates. ', 'normal'),
        (str(fail_num_dict['Gender']), 'bold'),
        (' sample(s) with discordant gender information are removed.', 'normal')
    ]
    
    ## content
    pdf.content(text_list)
    pdf.ln(10)

    # Missingness
    ## title
    pdf.chapter_title("4.2", "Individual level QC --- Missingness")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    ## text and format
    text_list = [
        ('Individual level Missingness is the number of SNPs that are missing for a specific individual. ', 'normal'),
        ('High levels of missingness indicates a technical problem or poor DNA quality within the data. ', 'normal'),
        ('The default threshold of missingness is set at ', 'normal'),
        (str(geno_miss), 'bold'),
        ('%. In this dataset, ', 'normal'),
        (str(fail_num_dict['Missingness']), 'bold'),
        (' individual(s) with high rates of genotype missingness were excluded.', 'normal')
    ]

    ## content
    pdf.content(text_list)
    pdf.ln(10)

    # Heterozygosity
    ## title
    pdf.chapter_title("4.3", "Individual level QC --- Heterozygosity")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    ## text and format
    text_list = [
        ('Heterozygosity is the carrying of two different alleles of a specific SNP. ', 'normal'),
        ('The heterozygosity rate of an individual is the proportion of heterozygous genotypes. ', 'normal'),
        ('High levels of heterozygosity within an individual might be an indication of low sample quality, ', 'normal'),
        ('whereas low levels of heterozygosity may be due to inbreeding. ', 'normal'),
        ('In this step, ', 'normal'),
        (str(fail_num_dict['Heterozygosity']), 'bold'),
        (' individuals who exceed ± ', 'normal'),
        (str(heter), 'bold'),
        (" SD from the samples' heterozygosity rate mean are excluded. ", 'normal')
    ]

    ## content
    pdf.content(text_list)
    pdf.ln(10)

    ## text and format
    text_list = [
        ('Heterozygosity and missingness are plotted on the next page based on all GWAS samples analysed in this study. ', 'normal'),
        ('The red horizontal lines indicate {0} standard deviations from the mean, '.format(heter), 'normal'),
        (' and the red vertical line indicates the missingness threshold of {0}'.format(geno_miss), 'normal')
    ]
    
    # content
    pdf.content(text_list)
    pdf.ln(10)
    
    # IBD
    ## titile
    pdf.chapter_title("4.4", "Individual level QC --- Identity By Descent (IBD) Test")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    ## text and format
    text_list = [
        ('The identity by descent (IBD) represents the degree of recent shared ancestry for a pair of individuals using a default cutoff set at ', 'normal'),
        (str(ibd), 'bold'),
        ('. IBD mapping relies on linkage disequilibrium (LD), which is used to keep a subset of SNPs that are seemingly uncorrelated with each other. ', 'normal'),
        ('In each LD region of ', 'normal'),
        (str(prune_window), 'bold'),
        (' variant counts, pairs of variants in the current window with a squared correlation greater than ', 'normal'),
        (str(prune_threshold), 'bold'),
        (' are removed. For each pair of individuals with high PI_HAT, the individual listed in the first column are excluded. Using this criteria, ', 'normal'),
        (str(fail_num_dict['IBD']), 'bold'),
        (' pairs of individuals failed the IBD test.', 'normal')
    ]

    ## content
    pdf.content(text_list)
    pdf.ln(10)

    # plot of het-vs-imiss
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(logo_file)
    pdf.image(args.fig_het_vs_imiss, x=10, w=190, h=180)

    # Population
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(logo_file)

    ## titile
    pdf.chapter_title("4.5", "Individual level QC --- Population Stratification")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    ## text and format
    text_list = [
        ('Populaiton stratification is the presence of a systematic difference in allele frequencies between subpopulaitons. ', 'normal'),
        (str(fail_num_dict['Population']), 'bold'),
        (' samples with different ethnic backgroud are removed. ', 'normal'),
        ('The following plot depicts the distribution of the data on the first two dimension computed by principal component analysis (PCA).', 'normal')
    ]

    ## content
    pdf.content(text_list)
    pdf.ln(10)

    ## figure
    pdf.image(args.fig_population, x=10, w=190, h=180)


    ################## SNPs QC #################################
    Title = 'SNPs level QC'
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(logo_file)

    # MAF
    ## title
    pdf.chapter_title("5.1", "SNPs level QC --- Minor Allele Frequency (MAF)")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    ## text and format
    text_list = [
        ('MAF is the frequency of the least often occurring allele at a specific location. ', 'normal'),
        ('Most studies are underpowered to detect associations with SNPs with a low MAF and therefore exclude these SNPs. ', 'normal'),
        ('In this study, the default MAF threshold is set at ', 'normal'),
        (str(maf), 'bold'),
        (', and ', 'normal'),
        (str(snp_fail_maf_num), 'bold'),
        (' SNPs with MAF lower than {0} are removed. '.format(maf), 'normal'),
        ('Plot shown below illustrates the histogram of MAF distribution, and the region with colored in red indicates the number of removed SNPs.', 'normal')
    ]

    ## content
    pdf.content(text_list)
    pdf.ln(10)

    ## figure
    pdf.image(args.fig_maf, x=10, w=190, h=180)

    # geno missingness
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(logo_file)

    ## title
    pdf.chapter_title("5.2", "SNPs level QC --- Missingness")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    ## text and format
    text_list = [
        ('SNP level missingness is the number of individuals in the sample for whom information on a specific SNP is missing. ', 'normal'),
        ('SNPs with a high level of missingness can potentially lead to bias. ', 'normal'),
        ('The threshold for SNP level missingness is set at ', 'normal'),
        (str(geno_miss), 'bold'),
        (', and ', 'normal'),
        (str(snp_fail_miss_num), 'bold'),
        (' SNPs are filtered out.', 'normal')
    ]

    ## content
    pdf.content(text_list)
    pdf.ln(10)

    ## figure
    pdf.image(args.fig_lmiss, x=10, w=190, h=180)

    # Hardy Weinberg Equilibrium
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(logo_file)

    ## title
    pdf.chapter_title("5.3", "SNPs level QC --- Hardy Weinberg Equilibrium (HWE)")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    ## text and format
    text_list = [
        ('HWE concerns the relationship between allele and genotype frequencies. ', 'normal'),
        ('Violation of the HWE can be attributed to inbreeding, contamination, or other forms of non-random mating. ', 'normal'),
        ('The threshold for SNP level missingness in this GWAS study is set at ', 'normal'),
        (str(hwe), 'bold'),
        (', resulting in ', 'normal'),
        (str(snp_fail_hwe_num), 'bold'),
        (' SNPs that violate the HWE law.', 'normal'),
        ('The following is the plot of the empirical cumulative distribution function (CDF), and the red region represents the percentage of removed SNPs.', 'normal')
    ]

    ## content
    pdf.content(text_list)
    pdf.ln(10)

    ## figure
    pdf.image(args.fig_hwe, x=10, w=190, h=180)


    ################## Report #################################
    pdf.output('{0}_qc_report.pdf'.format(basename), 'F')


if __name__ == "__main__":
    main()
