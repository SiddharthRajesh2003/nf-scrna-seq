#!/bin/bash

# Download and organize FASTQ files by sample
# Files are renamed to Cell Ranger format
# Naming scheme: sequential
# Uses SRA Toolkit (prefetch + fasterq-dump)

set -e

# Check if required tools are installed
command -v prefetch >/dev/null 2>&1 || { echo 'prefetch not found. Install SRA Toolkit.'; exit 1; }
command -v fasterq-dump >/dev/null 2>&1 || { echo 'fasterq-dump not found. Install SRA Toolkit.'; exit 1; }
command -v gzip >/dev/null 2>&1 || { echo 'gzip not found. Install gzip for compression.'; exit 1; }

# Sample Mapping:
# GSM6616991 -> NoOSA_001
# GSM6616992 -> NoOSA_002
# GSM6616993 -> NoOSA_003
# GSM6616994 -> NoOSA_004
# GSM6616995 -> NoOSA_005
# GSM6616996 -> NoOSA_006
# GSM6616997 -> NoOSA_007
# GSM6616998 -> NoOSA_008
# GSM6616999 -> NoOSA_009
# GSM6617000 -> NoOSA_010
# GSM6617001 -> NoOSA_011
# GSM6617002 -> OSA_001
# GSM6617003 -> OSA_002
# GSM6617004 -> OSA_003
# GSM6617005 -> OSA_004
# GSM6617006 -> OSA_005
# GSM6617007 -> OSA_006
# GSM6617008 -> OSA_007
# GSM6617009 -> OSA_008
# GSM6617010 -> OSA_009
# GSM6617011 -> OSA_010
# GSM6617012 -> OSA_011


#============================================================
# Sample: NoOSA_001 (Original: GSM6616991)
# Status: NoOSA, Sex: male, Age: 9
#============================================================

echo 'Processing NoOSA_001...'
mkdir -p fastq/NoOSA_001
cd fastq/NoOSA_001

# Lane 1/4: SRR21817835
echo '  Downloading SRR21817835 (Lane 1)...'
prefetch SRR21817835 || { echo 'Failed to prefetch SRR21817835'; exit 1; }
fasterq-dump SRR21817835 --split-files --threads 8 || { echo 'Failed to dump SRR21817835'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817835_1.fastq
gzip SRR21817835_2.fastq

# Rename to Cell Ranger format
mv SRR21817835_1.fastq.gz NoOSA_001_S1_L001_R1_001.fastq.gz
mv SRR21817835_2.fastq.gz NoOSA_001_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817835

# Lane 2/4: SRR21817836
echo '  Downloading SRR21817836 (Lane 2)...'
prefetch SRR21817836 || { echo 'Failed to prefetch SRR21817836'; exit 1; }
fasterq-dump SRR21817836 --split-files --threads 8 || { echo 'Failed to dump SRR21817836'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817836_1.fastq
gzip SRR21817836_2.fastq

# Rename to Cell Ranger format
mv SRR21817836_1.fastq.gz NoOSA_001_S1_L002_R1_001.fastq.gz
mv SRR21817836_2.fastq.gz NoOSA_001_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817836

# Lane 3/4: SRR21817837
echo '  Downloading SRR21817837 (Lane 3)...'
prefetch SRR21817837 || { echo 'Failed to prefetch SRR21817837'; exit 1; }
fasterq-dump SRR21817837 --split-files --threads 8 || { echo 'Failed to dump SRR21817837'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817837_1.fastq
gzip SRR21817837_2.fastq

# Rename to Cell Ranger format
mv SRR21817837_1.fastq.gz NoOSA_001_S1_L003_R1_001.fastq.gz
mv SRR21817837_2.fastq.gz NoOSA_001_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817837

# Lane 4/4: SRR21817838
echo '  Downloading SRR21817838 (Lane 4)...'
prefetch SRR21817838 || { echo 'Failed to prefetch SRR21817838'; exit 1; }
fasterq-dump SRR21817838 --split-files --threads 8 || { echo 'Failed to dump SRR21817838'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817838_1.fastq
gzip SRR21817838_2.fastq

# Rename to Cell Ranger format
mv SRR21817838_1.fastq.gz NoOSA_001_S1_L004_R1_001.fastq.gz
mv SRR21817838_2.fastq.gz NoOSA_001_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817838

cd ../..
echo '✓ Completed NoOSA_001'
echo ''


#============================================================
# Sample: NoOSA_002 (Original: GSM6616992)
# Status: NoOSA, Sex: male, Age: 4
#============================================================

echo 'Processing NoOSA_002...'
mkdir -p fastq/NoOSA_002
cd fastq/NoOSA_002

# Lane 1/4: SRR21817831
echo '  Downloading SRR21817831 (Lane 1)...'
prefetch SRR21817831 || { echo 'Failed to prefetch SRR21817831'; exit 1; }
fasterq-dump SRR21817831 --split-files --threads 8 || { echo 'Failed to dump SRR21817831'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817831_1.fastq
gzip SRR21817831_2.fastq

# Rename to Cell Ranger format
mv SRR21817831_1.fastq.gz NoOSA_002_S1_L001_R1_001.fastq.gz
mv SRR21817831_2.fastq.gz NoOSA_002_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817831

# Lane 2/4: SRR21817832
echo '  Downloading SRR21817832 (Lane 2)...'
prefetch SRR21817832 || { echo 'Failed to prefetch SRR21817832'; exit 1; }
fasterq-dump SRR21817832 --split-files --threads 8 || { echo 'Failed to dump SRR21817832'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817832_1.fastq
gzip SRR21817832_2.fastq

# Rename to Cell Ranger format
mv SRR21817832_1.fastq.gz NoOSA_002_S1_L002_R1_001.fastq.gz
mv SRR21817832_2.fastq.gz NoOSA_002_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817832

# Lane 3/4: SRR21817833
echo '  Downloading SRR21817833 (Lane 3)...'
prefetch SRR21817833 || { echo 'Failed to prefetch SRR21817833'; exit 1; }
fasterq-dump SRR21817833 --split-files --threads 8 || { echo 'Failed to dump SRR21817833'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817833_1.fastq
gzip SRR21817833_2.fastq

# Rename to Cell Ranger format
mv SRR21817833_1.fastq.gz NoOSA_002_S1_L003_R1_001.fastq.gz
mv SRR21817833_2.fastq.gz NoOSA_002_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817833

# Lane 4/4: SRR21817834
echo '  Downloading SRR21817834 (Lane 4)...'
prefetch SRR21817834 || { echo 'Failed to prefetch SRR21817834'; exit 1; }
fasterq-dump SRR21817834 --split-files --threads 8 || { echo 'Failed to dump SRR21817834'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817834_1.fastq
gzip SRR21817834_2.fastq

# Rename to Cell Ranger format
mv SRR21817834_1.fastq.gz NoOSA_002_S1_L004_R1_001.fastq.gz
mv SRR21817834_2.fastq.gz NoOSA_002_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817834

cd ../..
echo '✓ Completed NoOSA_002'
echo ''


#============================================================
# Sample: NoOSA_003 (Original: GSM6616993)
# Status: NoOSA, Sex: female, Age: 3
#============================================================

echo 'Processing NoOSA_003...'
mkdir -p fastq/NoOSA_003
cd fastq/NoOSA_003

# Lane 1/4: SRR21817827
echo '  Downloading SRR21817827 (Lane 1)...'
prefetch SRR21817827 || { echo 'Failed to prefetch SRR21817827'; exit 1; }
fasterq-dump SRR21817827 --split-files --threads 8 || { echo 'Failed to dump SRR21817827'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817827_1.fastq
gzip SRR21817827_2.fastq

# Rename to Cell Ranger format
mv SRR21817827_1.fastq.gz NoOSA_003_S1_L001_R1_001.fastq.gz
mv SRR21817827_2.fastq.gz NoOSA_003_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817827

# Lane 2/4: SRR21817828
echo '  Downloading SRR21817828 (Lane 2)...'
prefetch SRR21817828 || { echo 'Failed to prefetch SRR21817828'; exit 1; }
fasterq-dump SRR21817828 --split-files --threads 8 || { echo 'Failed to dump SRR21817828'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817828_1.fastq
gzip SRR21817828_2.fastq

# Rename to Cell Ranger format
mv SRR21817828_1.fastq.gz NoOSA_003_S1_L002_R1_001.fastq.gz
mv SRR21817828_2.fastq.gz NoOSA_003_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817828

# Lane 3/4: SRR21817829
echo '  Downloading SRR21817829 (Lane 3)...'
prefetch SRR21817829 || { echo 'Failed to prefetch SRR21817829'; exit 1; }
fasterq-dump SRR21817829 --split-files --threads 8 || { echo 'Failed to dump SRR21817829'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817829_1.fastq
gzip SRR21817829_2.fastq

# Rename to Cell Ranger format
mv SRR21817829_1.fastq.gz NoOSA_003_S1_L003_R1_001.fastq.gz
mv SRR21817829_2.fastq.gz NoOSA_003_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817829

# Lane 4/4: SRR21817830
echo '  Downloading SRR21817830 (Lane 4)...'
prefetch SRR21817830 || { echo 'Failed to prefetch SRR21817830'; exit 1; }
fasterq-dump SRR21817830 --split-files --threads 8 || { echo 'Failed to dump SRR21817830'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817830_1.fastq
gzip SRR21817830_2.fastq

# Rename to Cell Ranger format
mv SRR21817830_1.fastq.gz NoOSA_003_S1_L004_R1_001.fastq.gz
mv SRR21817830_2.fastq.gz NoOSA_003_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817830

cd ../..
echo '✓ Completed NoOSA_003'
echo ''


#============================================================
# Sample: NoOSA_004 (Original: GSM6616994)
# Status: NoOSA, Sex: male, Age: 8
#============================================================

echo 'Processing NoOSA_004...'
mkdir -p fastq/NoOSA_004
cd fastq/NoOSA_004

# Lane 1/4: SRR21817823
echo '  Downloading SRR21817823 (Lane 1)...'
prefetch SRR21817823 || { echo 'Failed to prefetch SRR21817823'; exit 1; }
fasterq-dump SRR21817823 --split-files --threads 8 || { echo 'Failed to dump SRR21817823'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817823_1.fastq
gzip SRR21817823_2.fastq

# Rename to Cell Ranger format
mv SRR21817823_1.fastq.gz NoOSA_004_S1_L001_R1_001.fastq.gz
mv SRR21817823_2.fastq.gz NoOSA_004_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817823

# Lane 2/4: SRR21817824
echo '  Downloading SRR21817824 (Lane 2)...'
prefetch SRR21817824 || { echo 'Failed to prefetch SRR21817824'; exit 1; }
fasterq-dump SRR21817824 --split-files --threads 8 || { echo 'Failed to dump SRR21817824'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817824_1.fastq
gzip SRR21817824_2.fastq

# Rename to Cell Ranger format
mv SRR21817824_1.fastq.gz NoOSA_004_S1_L002_R1_001.fastq.gz
mv SRR21817824_2.fastq.gz NoOSA_004_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817824

# Lane 3/4: SRR21817825
echo '  Downloading SRR21817825 (Lane 3)...'
prefetch SRR21817825 || { echo 'Failed to prefetch SRR21817825'; exit 1; }
fasterq-dump SRR21817825 --split-files --threads 8 || { echo 'Failed to dump SRR21817825'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817825_1.fastq
gzip SRR21817825_2.fastq

# Rename to Cell Ranger format
mv SRR21817825_1.fastq.gz NoOSA_004_S1_L003_R1_001.fastq.gz
mv SRR21817825_2.fastq.gz NoOSA_004_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817825

# Lane 4/4: SRR21817826
echo '  Downloading SRR21817826 (Lane 4)...'
prefetch SRR21817826 || { echo 'Failed to prefetch SRR21817826'; exit 1; }
fasterq-dump SRR21817826 --split-files --threads 8 || { echo 'Failed to dump SRR21817826'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817826_1.fastq
gzip SRR21817826_2.fastq

# Rename to Cell Ranger format
mv SRR21817826_1.fastq.gz NoOSA_004_S1_L004_R1_001.fastq.gz
mv SRR21817826_2.fastq.gz NoOSA_004_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817826

cd ../..
echo '✓ Completed NoOSA_004'
echo ''


#============================================================
# Sample: NoOSA_005 (Original: GSM6616995)
# Status: NoOSA, Sex: male, Age: 5
#============================================================

echo 'Processing NoOSA_005...'
mkdir -p fastq/NoOSA_005
cd fastq/NoOSA_005

# Lane 1/4: SRR21817819
echo '  Downloading SRR21817819 (Lane 1)...'
prefetch SRR21817819 || { echo 'Failed to prefetch SRR21817819'; exit 1; }
fasterq-dump SRR21817819 --split-files --threads 8 || { echo 'Failed to dump SRR21817819'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817819_1.fastq
gzip SRR21817819_2.fastq

# Rename to Cell Ranger format
mv SRR21817819_1.fastq.gz NoOSA_005_S1_L001_R1_001.fastq.gz
mv SRR21817819_2.fastq.gz NoOSA_005_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817819

# Lane 2/4: SRR21817820
echo '  Downloading SRR21817820 (Lane 2)...'
prefetch SRR21817820 || { echo 'Failed to prefetch SRR21817820'; exit 1; }
fasterq-dump SRR21817820 --split-files --threads 8 || { echo 'Failed to dump SRR21817820'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817820_1.fastq
gzip SRR21817820_2.fastq

# Rename to Cell Ranger format
mv SRR21817820_1.fastq.gz NoOSA_005_S1_L002_R1_001.fastq.gz
mv SRR21817820_2.fastq.gz NoOSA_005_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817820

# Lane 3/4: SRR21817821
echo '  Downloading SRR21817821 (Lane 3)...'
prefetch SRR21817821 || { echo 'Failed to prefetch SRR21817821'; exit 1; }
fasterq-dump SRR21817821 --split-files --threads 8 || { echo 'Failed to dump SRR21817821'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817821_1.fastq
gzip SRR21817821_2.fastq

# Rename to Cell Ranger format
mv SRR21817821_1.fastq.gz NoOSA_005_S1_L003_R1_001.fastq.gz
mv SRR21817821_2.fastq.gz NoOSA_005_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817821

# Lane 4/4: SRR21817822
echo '  Downloading SRR21817822 (Lane 4)...'
prefetch SRR21817822 || { echo 'Failed to prefetch SRR21817822'; exit 1; }
fasterq-dump SRR21817822 --split-files --threads 8 || { echo 'Failed to dump SRR21817822'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817822_1.fastq
gzip SRR21817822_2.fastq

# Rename to Cell Ranger format
mv SRR21817822_1.fastq.gz NoOSA_005_S1_L004_R1_001.fastq.gz
mv SRR21817822_2.fastq.gz NoOSA_005_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817822

cd ../..
echo '✓ Completed NoOSA_005'
echo ''


#============================================================
# Sample: NoOSA_006 (Original: GSM6616996)
# Status: NoOSA, Sex: female, Age: 5
#============================================================

echo 'Processing NoOSA_006...'
mkdir -p fastq/NoOSA_006
cd fastq/NoOSA_006

# Lane 1/4: SRR21817815
echo '  Downloading SRR21817815 (Lane 1)...'
prefetch SRR21817815 || { echo 'Failed to prefetch SRR21817815'; exit 1; }
fasterq-dump SRR21817815 --split-files --threads 8 || { echo 'Failed to dump SRR21817815'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817815_1.fastq
gzip SRR21817815_2.fastq

# Rename to Cell Ranger format
mv SRR21817815_1.fastq.gz NoOSA_006_S1_L001_R1_001.fastq.gz
mv SRR21817815_2.fastq.gz NoOSA_006_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817815

# Lane 2/4: SRR21817816
echo '  Downloading SRR21817816 (Lane 2)...'
prefetch SRR21817816 || { echo 'Failed to prefetch SRR21817816'; exit 1; }
fasterq-dump SRR21817816 --split-files --threads 8 || { echo 'Failed to dump SRR21817816'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817816_1.fastq
gzip SRR21817816_2.fastq

# Rename to Cell Ranger format
mv SRR21817816_1.fastq.gz NoOSA_006_S1_L002_R1_001.fastq.gz
mv SRR21817816_2.fastq.gz NoOSA_006_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817816

# Lane 3/4: SRR21817817
echo '  Downloading SRR21817817 (Lane 3)...'
prefetch SRR21817817 || { echo 'Failed to prefetch SRR21817817'; exit 1; }
fasterq-dump SRR21817817 --split-files --threads 8 || { echo 'Failed to dump SRR21817817'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817817_1.fastq
gzip SRR21817817_2.fastq

# Rename to Cell Ranger format
mv SRR21817817_1.fastq.gz NoOSA_006_S1_L003_R1_001.fastq.gz
mv SRR21817817_2.fastq.gz NoOSA_006_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817817

# Lane 4/4: SRR21817818
echo '  Downloading SRR21817818 (Lane 4)...'
prefetch SRR21817818 || { echo 'Failed to prefetch SRR21817818'; exit 1; }
fasterq-dump SRR21817818 --split-files --threads 8 || { echo 'Failed to dump SRR21817818'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817818_1.fastq
gzip SRR21817818_2.fastq

# Rename to Cell Ranger format
mv SRR21817818_1.fastq.gz NoOSA_006_S1_L004_R1_001.fastq.gz
mv SRR21817818_2.fastq.gz NoOSA_006_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817818

cd ../..
echo '✓ Completed NoOSA_006'
echo ''


#============================================================
# Sample: NoOSA_007 (Original: GSM6616997)
# Status: NoOSA, Sex: female, Age: 8
#============================================================

echo 'Processing NoOSA_007...'
mkdir -p fastq/NoOSA_007
cd fastq/NoOSA_007

# Lane 1/4: SRR21817811
echo '  Downloading SRR21817811 (Lane 1)...'
prefetch SRR21817811 || { echo 'Failed to prefetch SRR21817811'; exit 1; }
fasterq-dump SRR21817811 --split-files --threads 8 || { echo 'Failed to dump SRR21817811'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817811_1.fastq
gzip SRR21817811_2.fastq

# Rename to Cell Ranger format
mv SRR21817811_1.fastq.gz NoOSA_007_S1_L001_R1_001.fastq.gz
mv SRR21817811_2.fastq.gz NoOSA_007_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817811

# Lane 2/4: SRR21817812
echo '  Downloading SRR21817812 (Lane 2)...'
prefetch SRR21817812 || { echo 'Failed to prefetch SRR21817812'; exit 1; }
fasterq-dump SRR21817812 --split-files --threads 8 || { echo 'Failed to dump SRR21817812'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817812_1.fastq
gzip SRR21817812_2.fastq

# Rename to Cell Ranger format
mv SRR21817812_1.fastq.gz NoOSA_007_S1_L002_R1_001.fastq.gz
mv SRR21817812_2.fastq.gz NoOSA_007_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817812

# Lane 3/4: SRR21817813
echo '  Downloading SRR21817813 (Lane 3)...'
prefetch SRR21817813 || { echo 'Failed to prefetch SRR21817813'; exit 1; }
fasterq-dump SRR21817813 --split-files --threads 8 || { echo 'Failed to dump SRR21817813'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817813_1.fastq
gzip SRR21817813_2.fastq

# Rename to Cell Ranger format
mv SRR21817813_1.fastq.gz NoOSA_007_S1_L003_R1_001.fastq.gz
mv SRR21817813_2.fastq.gz NoOSA_007_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817813

# Lane 4/4: SRR21817814
echo '  Downloading SRR21817814 (Lane 4)...'
prefetch SRR21817814 || { echo 'Failed to prefetch SRR21817814'; exit 1; }
fasterq-dump SRR21817814 --split-files --threads 8 || { echo 'Failed to dump SRR21817814'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817814_1.fastq
gzip SRR21817814_2.fastq

# Rename to Cell Ranger format
mv SRR21817814_1.fastq.gz NoOSA_007_S1_L004_R1_001.fastq.gz
mv SRR21817814_2.fastq.gz NoOSA_007_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817814

cd ../..
echo '✓ Completed NoOSA_007'
echo ''


#============================================================
# Sample: NoOSA_008 (Original: GSM6616998)
# Status: NoOSA, Sex: female, Age: 5
#============================================================

echo 'Processing NoOSA_008...'
mkdir -p fastq/NoOSA_008
cd fastq/NoOSA_008

# Lane 1/4: SRR21817803
echo '  Downloading SRR21817803 (Lane 1)...'
prefetch SRR21817803 || { echo 'Failed to prefetch SRR21817803'; exit 1; }
fasterq-dump SRR21817803 --split-files --threads 8 || { echo 'Failed to dump SRR21817803'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817803_1.fastq
gzip SRR21817803_2.fastq

# Rename to Cell Ranger format
mv SRR21817803_1.fastq.gz NoOSA_008_S1_L001_R1_001.fastq.gz
mv SRR21817803_2.fastq.gz NoOSA_008_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817803

# Lane 2/4: SRR21817804
echo '  Downloading SRR21817804 (Lane 2)...'
prefetch SRR21817804 || { echo 'Failed to prefetch SRR21817804'; exit 1; }
fasterq-dump SRR21817804 --split-files --threads 8 || { echo 'Failed to dump SRR21817804'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817804_1.fastq
gzip SRR21817804_2.fastq

# Rename to Cell Ranger format
mv SRR21817804_1.fastq.gz NoOSA_008_S1_L002_R1_001.fastq.gz
mv SRR21817804_2.fastq.gz NoOSA_008_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817804

# Lane 3/4: SRR21817805
echo '  Downloading SRR21817805 (Lane 3)...'
prefetch SRR21817805 || { echo 'Failed to prefetch SRR21817805'; exit 1; }
fasterq-dump SRR21817805 --split-files --threads 8 || { echo 'Failed to dump SRR21817805'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817805_1.fastq
gzip SRR21817805_2.fastq

# Rename to Cell Ranger format
mv SRR21817805_1.fastq.gz NoOSA_008_S1_L003_R1_001.fastq.gz
mv SRR21817805_2.fastq.gz NoOSA_008_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817805

# Lane 4/4: SRR21817806
echo '  Downloading SRR21817806 (Lane 4)...'
prefetch SRR21817806 || { echo 'Failed to prefetch SRR21817806'; exit 1; }
fasterq-dump SRR21817806 --split-files --threads 8 || { echo 'Failed to dump SRR21817806'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817806_1.fastq
gzip SRR21817806_2.fastq

# Rename to Cell Ranger format
mv SRR21817806_1.fastq.gz NoOSA_008_S1_L004_R1_001.fastq.gz
mv SRR21817806_2.fastq.gz NoOSA_008_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817806

cd ../..
echo '✓ Completed NoOSA_008'
echo ''


#============================================================
# Sample: NoOSA_009 (Original: GSM6616999)
# Status: NoOSA, Sex: male, Age: 5
#============================================================

echo 'Processing NoOSA_009...'
mkdir -p fastq/NoOSA_009
cd fastq/NoOSA_009

# Lane 1/4: SRR21817795
echo '  Downloading SRR21817795 (Lane 1)...'
prefetch SRR21817795 || { echo 'Failed to prefetch SRR21817795'; exit 1; }
fasterq-dump SRR21817795 --split-files --threads 8 || { echo 'Failed to dump SRR21817795'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817795_1.fastq
gzip SRR21817795_2.fastq

# Rename to Cell Ranger format
mv SRR21817795_1.fastq.gz NoOSA_009_S1_L001_R1_001.fastq.gz
mv SRR21817795_2.fastq.gz NoOSA_009_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817795

# Lane 2/4: SRR21817796
echo '  Downloading SRR21817796 (Lane 2)...'
prefetch SRR21817796 || { echo 'Failed to prefetch SRR21817796'; exit 1; }
fasterq-dump SRR21817796 --split-files --threads 8 || { echo 'Failed to dump SRR21817796'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817796_1.fastq
gzip SRR21817796_2.fastq

# Rename to Cell Ranger format
mv SRR21817796_1.fastq.gz NoOSA_009_S1_L002_R1_001.fastq.gz
mv SRR21817796_2.fastq.gz NoOSA_009_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817796

# Lane 3/4: SRR21817797
echo '  Downloading SRR21817797 (Lane 3)...'
prefetch SRR21817797 || { echo 'Failed to prefetch SRR21817797'; exit 1; }
fasterq-dump SRR21817797 --split-files --threads 8 || { echo 'Failed to dump SRR21817797'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817797_1.fastq
gzip SRR21817797_2.fastq

# Rename to Cell Ranger format
mv SRR21817797_1.fastq.gz NoOSA_009_S1_L003_R1_001.fastq.gz
mv SRR21817797_2.fastq.gz NoOSA_009_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817797

# Lane 4/4: SRR21817798
echo '  Downloading SRR21817798 (Lane 4)...'
prefetch SRR21817798 || { echo 'Failed to prefetch SRR21817798'; exit 1; }
fasterq-dump SRR21817798 --split-files --threads 8 || { echo 'Failed to dump SRR21817798'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817798_1.fastq
gzip SRR21817798_2.fastq

# Rename to Cell Ranger format
mv SRR21817798_1.fastq.gz NoOSA_009_S1_L004_R1_001.fastq.gz
mv SRR21817798_2.fastq.gz NoOSA_009_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817798

cd ../..
echo '✓ Completed NoOSA_009'
echo ''


#============================================================
# Sample: NoOSA_010 (Original: GSM6617000)
# Status: NoOSA, Sex: male, Age: 7
#============================================================

echo 'Processing NoOSA_010...'
mkdir -p fastq/NoOSA_010
cd fastq/NoOSA_010

# Lane 1/4: SRR21817791
echo '  Downloading SRR21817791 (Lane 1)...'
prefetch SRR21817791 || { echo 'Failed to prefetch SRR21817791'; exit 1; }
fasterq-dump SRR21817791 --split-files --threads 8 || { echo 'Failed to dump SRR21817791'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817791_1.fastq
gzip SRR21817791_2.fastq

# Rename to Cell Ranger format
mv SRR21817791_1.fastq.gz NoOSA_010_S1_L001_R1_001.fastq.gz
mv SRR21817791_2.fastq.gz NoOSA_010_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817791

# Lane 2/4: SRR21817792
echo '  Downloading SRR21817792 (Lane 2)...'
prefetch SRR21817792 || { echo 'Failed to prefetch SRR21817792'; exit 1; }
fasterq-dump SRR21817792 --split-files --threads 8 || { echo 'Failed to dump SRR21817792'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817792_1.fastq
gzip SRR21817792_2.fastq

# Rename to Cell Ranger format
mv SRR21817792_1.fastq.gz NoOSA_010_S1_L002_R1_001.fastq.gz
mv SRR21817792_2.fastq.gz NoOSA_010_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817792

# Lane 3/4: SRR21817793
echo '  Downloading SRR21817793 (Lane 3)...'
prefetch SRR21817793 || { echo 'Failed to prefetch SRR21817793'; exit 1; }
fasterq-dump SRR21817793 --split-files --threads 8 || { echo 'Failed to dump SRR21817793'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817793_1.fastq
gzip SRR21817793_2.fastq

# Rename to Cell Ranger format
mv SRR21817793_1.fastq.gz NoOSA_010_S1_L003_R1_001.fastq.gz
mv SRR21817793_2.fastq.gz NoOSA_010_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817793

# Lane 4/4: SRR21817794
echo '  Downloading SRR21817794 (Lane 4)...'
prefetch SRR21817794 || { echo 'Failed to prefetch SRR21817794'; exit 1; }
fasterq-dump SRR21817794 --split-files --threads 8 || { echo 'Failed to dump SRR21817794'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817794_1.fastq
gzip SRR21817794_2.fastq

# Rename to Cell Ranger format
mv SRR21817794_1.fastq.gz NoOSA_010_S1_L004_R1_001.fastq.gz
mv SRR21817794_2.fastq.gz NoOSA_010_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817794

cd ../..
echo '✓ Completed NoOSA_010'
echo ''


#============================================================
# Sample: NoOSA_011 (Original: GSM6617001)
# Status: NoOSA, Sex: male, Age: 4
#============================================================

echo 'Processing NoOSA_011...'
mkdir -p fastq/NoOSA_011
cd fastq/NoOSA_011

# Lane 1/4: SRR21817787
echo '  Downloading SRR21817787 (Lane 1)...'
prefetch SRR21817787 || { echo 'Failed to prefetch SRR21817787'; exit 1; }
fasterq-dump SRR21817787 --split-files --threads 8 || { echo 'Failed to dump SRR21817787'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817787_1.fastq
gzip SRR21817787_2.fastq

# Rename to Cell Ranger format
mv SRR21817787_1.fastq.gz NoOSA_011_S1_L001_R1_001.fastq.gz
mv SRR21817787_2.fastq.gz NoOSA_011_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817787

# Lane 2/4: SRR21817788
echo '  Downloading SRR21817788 (Lane 2)...'
prefetch SRR21817788 || { echo 'Failed to prefetch SRR21817788'; exit 1; }
fasterq-dump SRR21817788 --split-files --threads 8 || { echo 'Failed to dump SRR21817788'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817788_1.fastq
gzip SRR21817788_2.fastq

# Rename to Cell Ranger format
mv SRR21817788_1.fastq.gz NoOSA_011_S1_L002_R1_001.fastq.gz
mv SRR21817788_2.fastq.gz NoOSA_011_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817788

# Lane 3/4: SRR21817789
echo '  Downloading SRR21817789 (Lane 3)...'
prefetch SRR21817789 || { echo 'Failed to prefetch SRR21817789'; exit 1; }
fasterq-dump SRR21817789 --split-files --threads 8 || { echo 'Failed to dump SRR21817789'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817789_1.fastq
gzip SRR21817789_2.fastq

# Rename to Cell Ranger format
mv SRR21817789_1.fastq.gz NoOSA_011_S1_L003_R1_001.fastq.gz
mv SRR21817789_2.fastq.gz NoOSA_011_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817789

# Lane 4/4: SRR21817790
echo '  Downloading SRR21817790 (Lane 4)...'
prefetch SRR21817790 || { echo 'Failed to prefetch SRR21817790'; exit 1; }
fasterq-dump SRR21817790 --split-files --threads 8 || { echo 'Failed to dump SRR21817790'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817790_1.fastq
gzip SRR21817790_2.fastq

# Rename to Cell Ranger format
mv SRR21817790_1.fastq.gz NoOSA_011_S1_L004_R1_001.fastq.gz
mv SRR21817790_2.fastq.gz NoOSA_011_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817790

cd ../..
echo '✓ Completed NoOSA_011'
echo ''


#============================================================
# Sample: OSA_001 (Original: GSM6617002)
# Status: OSA, Sex: male, Age: 2
#============================================================

echo 'Processing OSA_001...'
mkdir -p fastq/OSA_001
cd fastq/OSA_001

# Lane 1/4: SRR21817783
echo '  Downloading SRR21817783 (Lane 1)...'
prefetch SRR21817783 || { echo 'Failed to prefetch SRR21817783'; exit 1; }
fasterq-dump SRR21817783 --split-files --threads 8 || { echo 'Failed to dump SRR21817783'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817783_1.fastq
gzip SRR21817783_2.fastq

# Rename to Cell Ranger format
mv SRR21817783_1.fastq.gz OSA_001_S1_L001_R1_001.fastq.gz
mv SRR21817783_2.fastq.gz OSA_001_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817783

# Lane 2/4: SRR21817784
echo '  Downloading SRR21817784 (Lane 2)...'
prefetch SRR21817784 || { echo 'Failed to prefetch SRR21817784'; exit 1; }
fasterq-dump SRR21817784 --split-files --threads 8 || { echo 'Failed to dump SRR21817784'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817784_1.fastq
gzip SRR21817784_2.fastq

# Rename to Cell Ranger format
mv SRR21817784_1.fastq.gz OSA_001_S1_L002_R1_001.fastq.gz
mv SRR21817784_2.fastq.gz OSA_001_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817784

# Lane 3/4: SRR21817785
echo '  Downloading SRR21817785 (Lane 3)...'
prefetch SRR21817785 || { echo 'Failed to prefetch SRR21817785'; exit 1; }
fasterq-dump SRR21817785 --split-files --threads 8 || { echo 'Failed to dump SRR21817785'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817785_1.fastq
gzip SRR21817785_2.fastq

# Rename to Cell Ranger format
mv SRR21817785_1.fastq.gz OSA_001_S1_L003_R1_001.fastq.gz
mv SRR21817785_2.fastq.gz OSA_001_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817785

# Lane 4/4: SRR21817786
echo '  Downloading SRR21817786 (Lane 4)...'
prefetch SRR21817786 || { echo 'Failed to prefetch SRR21817786'; exit 1; }
fasterq-dump SRR21817786 --split-files --threads 8 || { echo 'Failed to dump SRR21817786'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817786_1.fastq
gzip SRR21817786_2.fastq

# Rename to Cell Ranger format
mv SRR21817786_1.fastq.gz OSA_001_S1_L004_R1_001.fastq.gz
mv SRR21817786_2.fastq.gz OSA_001_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817786

cd ../..
echo '✓ Completed OSA_001'
echo ''


#============================================================
# Sample: OSA_002 (Original: GSM6617003)
# Status: OSA, Sex: female, Age: 12
#============================================================

echo 'Processing OSA_002...'
mkdir -p fastq/OSA_002
cd fastq/OSA_002

# Lane 1/4: SRR21817779
echo '  Downloading SRR21817779 (Lane 1)...'
prefetch SRR21817779 || { echo 'Failed to prefetch SRR21817779'; exit 1; }
fasterq-dump SRR21817779 --split-files --threads 8 || { echo 'Failed to dump SRR21817779'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817779_1.fastq
gzip SRR21817779_2.fastq

# Rename to Cell Ranger format
mv SRR21817779_1.fastq.gz OSA_002_S1_L001_R1_001.fastq.gz
mv SRR21817779_2.fastq.gz OSA_002_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817779

# Lane 2/4: SRR21817780
echo '  Downloading SRR21817780 (Lane 2)...'
prefetch SRR21817780 || { echo 'Failed to prefetch SRR21817780'; exit 1; }
fasterq-dump SRR21817780 --split-files --threads 8 || { echo 'Failed to dump SRR21817780'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817780_1.fastq
gzip SRR21817780_2.fastq

# Rename to Cell Ranger format
mv SRR21817780_1.fastq.gz OSA_002_S1_L002_R1_001.fastq.gz
mv SRR21817780_2.fastq.gz OSA_002_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817780

# Lane 3/4: SRR21817781
echo '  Downloading SRR21817781 (Lane 3)...'
prefetch SRR21817781 || { echo 'Failed to prefetch SRR21817781'; exit 1; }
fasterq-dump SRR21817781 --split-files --threads 8 || { echo 'Failed to dump SRR21817781'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817781_1.fastq
gzip SRR21817781_2.fastq

# Rename to Cell Ranger format
mv SRR21817781_1.fastq.gz OSA_002_S1_L003_R1_001.fastq.gz
mv SRR21817781_2.fastq.gz OSA_002_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817781

# Lane 4/4: SRR21817782
echo '  Downloading SRR21817782 (Lane 4)...'
prefetch SRR21817782 || { echo 'Failed to prefetch SRR21817782'; exit 1; }
fasterq-dump SRR21817782 --split-files --threads 8 || { echo 'Failed to dump SRR21817782'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817782_1.fastq
gzip SRR21817782_2.fastq

# Rename to Cell Ranger format
mv SRR21817782_1.fastq.gz OSA_002_S1_L004_R1_001.fastq.gz
mv SRR21817782_2.fastq.gz OSA_002_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817782

cd ../..
echo '✓ Completed OSA_002'
echo ''


#============================================================
# Sample: OSA_003 (Original: GSM6617004)
# Status: OSA, Sex: male, Age: 7
#============================================================

echo 'Processing OSA_003...'
mkdir -p fastq/OSA_003
cd fastq/OSA_003

# Lane 1/4: SRR21817775
echo '  Downloading SRR21817775 (Lane 1)...'
prefetch SRR21817775 || { echo 'Failed to prefetch SRR21817775'; exit 1; }
fasterq-dump SRR21817775 --split-files --threads 8 || { echo 'Failed to dump SRR21817775'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817775_1.fastq
gzip SRR21817775_2.fastq

# Rename to Cell Ranger format
mv SRR21817775_1.fastq.gz OSA_003_S1_L001_R1_001.fastq.gz
mv SRR21817775_2.fastq.gz OSA_003_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817775

# Lane 2/4: SRR21817776
echo '  Downloading SRR21817776 (Lane 2)...'
prefetch SRR21817776 || { echo 'Failed to prefetch SRR21817776'; exit 1; }
fasterq-dump SRR21817776 --split-files --threads 8 || { echo 'Failed to dump SRR21817776'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817776_1.fastq
gzip SRR21817776_2.fastq

# Rename to Cell Ranger format
mv SRR21817776_1.fastq.gz OSA_003_S1_L002_R1_001.fastq.gz
mv SRR21817776_2.fastq.gz OSA_003_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817776

# Lane 3/4: SRR21817777
echo '  Downloading SRR21817777 (Lane 3)...'
prefetch SRR21817777 || { echo 'Failed to prefetch SRR21817777'; exit 1; }
fasterq-dump SRR21817777 --split-files --threads 8 || { echo 'Failed to dump SRR21817777'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817777_1.fastq
gzip SRR21817777_2.fastq

# Rename to Cell Ranger format
mv SRR21817777_1.fastq.gz OSA_003_S1_L003_R1_001.fastq.gz
mv SRR21817777_2.fastq.gz OSA_003_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817777

# Lane 4/4: SRR21817778
echo '  Downloading SRR21817778 (Lane 4)...'
prefetch SRR21817778 || { echo 'Failed to prefetch SRR21817778'; exit 1; }
fasterq-dump SRR21817778 --split-files --threads 8 || { echo 'Failed to dump SRR21817778'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817778_1.fastq
gzip SRR21817778_2.fastq

# Rename to Cell Ranger format
mv SRR21817778_1.fastq.gz OSA_003_S1_L004_R1_001.fastq.gz
mv SRR21817778_2.fastq.gz OSA_003_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817778

cd ../..
echo '✓ Completed OSA_003'
echo ''


#============================================================
# Sample: OSA_004 (Original: GSM6617005)
# Status: OSA, Sex: female, Age: 4
#============================================================

echo 'Processing OSA_004...'
mkdir -p fastq/OSA_004
cd fastq/OSA_004

# Lane 1/4: SRR21817771
echo '  Downloading SRR21817771 (Lane 1)...'
prefetch SRR21817771 || { echo 'Failed to prefetch SRR21817771'; exit 1; }
fasterq-dump SRR21817771 --split-files --threads 8 || { echo 'Failed to dump SRR21817771'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817771_1.fastq
gzip SRR21817771_2.fastq

# Rename to Cell Ranger format
mv SRR21817771_1.fastq.gz OSA_004_S1_L001_R1_001.fastq.gz
mv SRR21817771_2.fastq.gz OSA_004_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817771

# Lane 2/4: SRR21817772
echo '  Downloading SRR21817772 (Lane 2)...'
prefetch SRR21817772 || { echo 'Failed to prefetch SRR21817772'; exit 1; }
fasterq-dump SRR21817772 --split-files --threads 8 || { echo 'Failed to dump SRR21817772'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817772_1.fastq
gzip SRR21817772_2.fastq

# Rename to Cell Ranger format
mv SRR21817772_1.fastq.gz OSA_004_S1_L002_R1_001.fastq.gz
mv SRR21817772_2.fastq.gz OSA_004_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817772

# Lane 3/4: SRR21817773
echo '  Downloading SRR21817773 (Lane 3)...'
prefetch SRR21817773 || { echo 'Failed to prefetch SRR21817773'; exit 1; }
fasterq-dump SRR21817773 --split-files --threads 8 || { echo 'Failed to dump SRR21817773'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817773_1.fastq
gzip SRR21817773_2.fastq

# Rename to Cell Ranger format
mv SRR21817773_1.fastq.gz OSA_004_S1_L003_R1_001.fastq.gz
mv SRR21817773_2.fastq.gz OSA_004_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817773

# Lane 4/4: SRR21817774
echo '  Downloading SRR21817774 (Lane 4)...'
prefetch SRR21817774 || { echo 'Failed to prefetch SRR21817774'; exit 1; }
fasterq-dump SRR21817774 --split-files --threads 8 || { echo 'Failed to dump SRR21817774'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817774_1.fastq
gzip SRR21817774_2.fastq

# Rename to Cell Ranger format
mv SRR21817774_1.fastq.gz OSA_004_S1_L004_R1_001.fastq.gz
mv SRR21817774_2.fastq.gz OSA_004_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817774

cd ../..
echo '✓ Completed OSA_004'
echo ''


#============================================================
# Sample: OSA_005 (Original: GSM6617006)
# Status: OSA, Sex: female, Age: 4
#============================================================

echo 'Processing OSA_005...'
mkdir -p fastq/OSA_005
cd fastq/OSA_005

# Lane 1/4: SRR21817767
echo '  Downloading SRR21817767 (Lane 1)...'
prefetch SRR21817767 || { echo 'Failed to prefetch SRR21817767'; exit 1; }
fasterq-dump SRR21817767 --split-files --threads 8 || { echo 'Failed to dump SRR21817767'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817767_1.fastq
gzip SRR21817767_2.fastq

# Rename to Cell Ranger format
mv SRR21817767_1.fastq.gz OSA_005_S1_L001_R1_001.fastq.gz
mv SRR21817767_2.fastq.gz OSA_005_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817767

# Lane 2/4: SRR21817768
echo '  Downloading SRR21817768 (Lane 2)...'
prefetch SRR21817768 || { echo 'Failed to prefetch SRR21817768'; exit 1; }
fasterq-dump SRR21817768 --split-files --threads 8 || { echo 'Failed to dump SRR21817768'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817768_1.fastq
gzip SRR21817768_2.fastq

# Rename to Cell Ranger format
mv SRR21817768_1.fastq.gz OSA_005_S1_L002_R1_001.fastq.gz
mv SRR21817768_2.fastq.gz OSA_005_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817768

# Lane 3/4: SRR21817769
echo '  Downloading SRR21817769 (Lane 3)...'
prefetch SRR21817769 || { echo 'Failed to prefetch SRR21817769'; exit 1; }
fasterq-dump SRR21817769 --split-files --threads 8 || { echo 'Failed to dump SRR21817769'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817769_1.fastq
gzip SRR21817769_2.fastq

# Rename to Cell Ranger format
mv SRR21817769_1.fastq.gz OSA_005_S1_L003_R1_001.fastq.gz
mv SRR21817769_2.fastq.gz OSA_005_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817769

# Lane 4/4: SRR21817770
echo '  Downloading SRR21817770 (Lane 4)...'
prefetch SRR21817770 || { echo 'Failed to prefetch SRR21817770'; exit 1; }
fasterq-dump SRR21817770 --split-files --threads 8 || { echo 'Failed to dump SRR21817770'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817770_1.fastq
gzip SRR21817770_2.fastq

# Rename to Cell Ranger format
mv SRR21817770_1.fastq.gz OSA_005_S1_L004_R1_001.fastq.gz
mv SRR21817770_2.fastq.gz OSA_005_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817770

cd ../..
echo '✓ Completed OSA_005'
echo ''


#============================================================
# Sample: OSA_006 (Original: GSM6617007)
# Status: OSA, Sex: male, Age: 9
#============================================================

echo 'Processing OSA_006...'
mkdir -p fastq/OSA_006
cd fastq/OSA_006

# Lane 1/4: SRR21817763
echo '  Downloading SRR21817763 (Lane 1)...'
prefetch SRR21817763 || { echo 'Failed to prefetch SRR21817763'; exit 1; }
fasterq-dump SRR21817763 --split-files --threads 8 || { echo 'Failed to dump SRR21817763'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817763_1.fastq
gzip SRR21817763_2.fastq

# Rename to Cell Ranger format
mv SRR21817763_1.fastq.gz OSA_006_S1_L001_R1_001.fastq.gz
mv SRR21817763_2.fastq.gz OSA_006_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817763

# Lane 2/4: SRR21817764
echo '  Downloading SRR21817764 (Lane 2)...'
prefetch SRR21817764 || { echo 'Failed to prefetch SRR21817764'; exit 1; }
fasterq-dump SRR21817764 --split-files --threads 8 || { echo 'Failed to dump SRR21817764'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817764_1.fastq
gzip SRR21817764_2.fastq

# Rename to Cell Ranger format
mv SRR21817764_1.fastq.gz OSA_006_S1_L002_R1_001.fastq.gz
mv SRR21817764_2.fastq.gz OSA_006_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817764

# Lane 3/4: SRR21817765
echo '  Downloading SRR21817765 (Lane 3)...'
prefetch SRR21817765 || { echo 'Failed to prefetch SRR21817765'; exit 1; }
fasterq-dump SRR21817765 --split-files --threads 8 || { echo 'Failed to dump SRR21817765'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817765_1.fastq
gzip SRR21817765_2.fastq

# Rename to Cell Ranger format
mv SRR21817765_1.fastq.gz OSA_006_S1_L003_R1_001.fastq.gz
mv SRR21817765_2.fastq.gz OSA_006_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817765

# Lane 4/4: SRR21817766
echo '  Downloading SRR21817766 (Lane 4)...'
prefetch SRR21817766 || { echo 'Failed to prefetch SRR21817766'; exit 1; }
fasterq-dump SRR21817766 --split-files --threads 8 || { echo 'Failed to dump SRR21817766'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817766_1.fastq
gzip SRR21817766_2.fastq

# Rename to Cell Ranger format
mv SRR21817766_1.fastq.gz OSA_006_S1_L004_R1_001.fastq.gz
mv SRR21817766_2.fastq.gz OSA_006_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817766

cd ../..
echo '✓ Completed OSA_006'
echo ''


#============================================================
# Sample: OSA_007 (Original: GSM6617008)
# Status: OSA, Sex: female, Age: 4
#============================================================

echo 'Processing OSA_007...'
mkdir -p fastq/OSA_007
cd fastq/OSA_007

# Lane 1/4: SRR21817759
echo '  Downloading SRR21817759 (Lane 1)...'
prefetch SRR21817759 || { echo 'Failed to prefetch SRR21817759'; exit 1; }
fasterq-dump SRR21817759 --split-files --threads 8 || { echo 'Failed to dump SRR21817759'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817759_1.fastq
gzip SRR21817759_2.fastq

# Rename to Cell Ranger format
mv SRR21817759_1.fastq.gz OSA_007_S1_L001_R1_001.fastq.gz
mv SRR21817759_2.fastq.gz OSA_007_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817759

# Lane 2/4: SRR21817760
echo '  Downloading SRR21817760 (Lane 2)...'
prefetch SRR21817760 || { echo 'Failed to prefetch SRR21817760'; exit 1; }
fasterq-dump SRR21817760 --split-files --threads 8 || { echo 'Failed to dump SRR21817760'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817760_1.fastq
gzip SRR21817760_2.fastq

# Rename to Cell Ranger format
mv SRR21817760_1.fastq.gz OSA_007_S1_L002_R1_001.fastq.gz
mv SRR21817760_2.fastq.gz OSA_007_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817760

# Lane 3/4: SRR21817761
echo '  Downloading SRR21817761 (Lane 3)...'
prefetch SRR21817761 || { echo 'Failed to prefetch SRR21817761'; exit 1; }
fasterq-dump SRR21817761 --split-files --threads 8 || { echo 'Failed to dump SRR21817761'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817761_1.fastq
gzip SRR21817761_2.fastq

# Rename to Cell Ranger format
mv SRR21817761_1.fastq.gz OSA_007_S1_L003_R1_001.fastq.gz
mv SRR21817761_2.fastq.gz OSA_007_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817761

# Lane 4/4: SRR21817762
echo '  Downloading SRR21817762 (Lane 4)...'
prefetch SRR21817762 || { echo 'Failed to prefetch SRR21817762'; exit 1; }
fasterq-dump SRR21817762 --split-files --threads 8 || { echo 'Failed to dump SRR21817762'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817762_1.fastq
gzip SRR21817762_2.fastq

# Rename to Cell Ranger format
mv SRR21817762_1.fastq.gz OSA_007_S1_L004_R1_001.fastq.gz
mv SRR21817762_2.fastq.gz OSA_007_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817762

cd ../..
echo '✓ Completed OSA_007'
echo ''


#============================================================
# Sample: OSA_008 (Original: GSM6617009)
# Status: OSA, Sex: female, Age: 4
#============================================================

echo 'Processing OSA_008...'
mkdir -p fastq/OSA_008
cd fastq/OSA_008

# Lane 1/4: SRR21817755
echo '  Downloading SRR21817755 (Lane 1)...'
prefetch SRR21817755 || { echo 'Failed to prefetch SRR21817755'; exit 1; }
fasterq-dump SRR21817755 --split-files --threads 8 || { echo 'Failed to dump SRR21817755'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817755_1.fastq
gzip SRR21817755_2.fastq

# Rename to Cell Ranger format
mv SRR21817755_1.fastq.gz OSA_008_S1_L001_R1_001.fastq.gz
mv SRR21817755_2.fastq.gz OSA_008_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817755

# Lane 2/4: SRR21817756
echo '  Downloading SRR21817756 (Lane 2)...'
prefetch SRR21817756 || { echo 'Failed to prefetch SRR21817756'; exit 1; }
fasterq-dump SRR21817756 --split-files --threads 8 || { echo 'Failed to dump SRR21817756'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817756_1.fastq
gzip SRR21817756_2.fastq

# Rename to Cell Ranger format
mv SRR21817756_1.fastq.gz OSA_008_S1_L002_R1_001.fastq.gz
mv SRR21817756_2.fastq.gz OSA_008_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817756

# Lane 3/4: SRR21817757
echo '  Downloading SRR21817757 (Lane 3)...'
prefetch SRR21817757 || { echo 'Failed to prefetch SRR21817757'; exit 1; }
fasterq-dump SRR21817757 --split-files --threads 8 || { echo 'Failed to dump SRR21817757'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817757_1.fastq
gzip SRR21817757_2.fastq

# Rename to Cell Ranger format
mv SRR21817757_1.fastq.gz OSA_008_S1_L003_R1_001.fastq.gz
mv SRR21817757_2.fastq.gz OSA_008_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817757

# Lane 4/4: SRR21817758
echo '  Downloading SRR21817758 (Lane 4)...'
prefetch SRR21817758 || { echo 'Failed to prefetch SRR21817758'; exit 1; }
fasterq-dump SRR21817758 --split-files --threads 8 || { echo 'Failed to dump SRR21817758'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817758_1.fastq
gzip SRR21817758_2.fastq

# Rename to Cell Ranger format
mv SRR21817758_1.fastq.gz OSA_008_S1_L004_R1_001.fastq.gz
mv SRR21817758_2.fastq.gz OSA_008_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817758

cd ../..
echo '✓ Completed OSA_008'
echo ''


#============================================================
# Sample: OSA_009 (Original: GSM6617010)
# Status: OSA, Sex: female, Age: 4
#============================================================

echo 'Processing OSA_009...'
mkdir -p fastq/OSA_009
cd fastq/OSA_009

# Lane 1/4: SRR21817751
echo '  Downloading SRR21817751 (Lane 1)...'
prefetch SRR21817751 || { echo 'Failed to prefetch SRR21817751'; exit 1; }
fasterq-dump SRR21817751 --split-files --threads 8 || { echo 'Failed to dump SRR21817751'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817751_1.fastq
gzip SRR21817751_2.fastq

# Rename to Cell Ranger format
mv SRR21817751_1.fastq.gz OSA_009_S1_L001_R1_001.fastq.gz
mv SRR21817751_2.fastq.gz OSA_009_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817751

# Lane 2/4: SRR21817752
echo '  Downloading SRR21817752 (Lane 2)...'
prefetch SRR21817752 || { echo 'Failed to prefetch SRR21817752'; exit 1; }
fasterq-dump SRR21817752 --split-files --threads 8 || { echo 'Failed to dump SRR21817752'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817752_1.fastq
gzip SRR21817752_2.fastq

# Rename to Cell Ranger format
mv SRR21817752_1.fastq.gz OSA_009_S1_L002_R1_001.fastq.gz
mv SRR21817752_2.fastq.gz OSA_009_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817752

# Lane 3/4: SRR21817753
echo '  Downloading SRR21817753 (Lane 3)...'
prefetch SRR21817753 || { echo 'Failed to prefetch SRR21817753'; exit 1; }
fasterq-dump SRR21817753 --split-files --threads 8 || { echo 'Failed to dump SRR21817753'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817753_1.fastq
gzip SRR21817753_2.fastq

# Rename to Cell Ranger format
mv SRR21817753_1.fastq.gz OSA_009_S1_L003_R1_001.fastq.gz
mv SRR21817753_2.fastq.gz OSA_009_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817753

# Lane 4/4: SRR21817754
echo '  Downloading SRR21817754 (Lane 4)...'
prefetch SRR21817754 || { echo 'Failed to prefetch SRR21817754'; exit 1; }
fasterq-dump SRR21817754 --split-files --threads 8 || { echo 'Failed to dump SRR21817754'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817754_1.fastq
gzip SRR21817754_2.fastq

# Rename to Cell Ranger format
mv SRR21817754_1.fastq.gz OSA_009_S1_L004_R1_001.fastq.gz
mv SRR21817754_2.fastq.gz OSA_009_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817754

cd ../..
echo '✓ Completed OSA_009'
echo ''


#============================================================
# Sample: OSA_010 (Original: GSM6617011)
# Status: OSA, Sex: male, Age: 16
#============================================================

echo 'Processing OSA_010...'
mkdir -p fastq/OSA_010
cd fastq/OSA_010

# Lane 1/4: SRR21817807
echo '  Downloading SRR21817807 (Lane 1)...'
prefetch SRR21817807 || { echo 'Failed to prefetch SRR21817807'; exit 1; }
fasterq-dump SRR21817807 --split-files --threads 8 || { echo 'Failed to dump SRR21817807'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817807_1.fastq
gzip SRR21817807_2.fastq

# Rename to Cell Ranger format
mv SRR21817807_1.fastq.gz OSA_010_S1_L001_R1_001.fastq.gz
mv SRR21817807_2.fastq.gz OSA_010_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817807

# Lane 2/4: SRR21817808
echo '  Downloading SRR21817808 (Lane 2)...'
prefetch SRR21817808 || { echo 'Failed to prefetch SRR21817808'; exit 1; }
fasterq-dump SRR21817808 --split-files --threads 8 || { echo 'Failed to dump SRR21817808'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817808_1.fastq
gzip SRR21817808_2.fastq

# Rename to Cell Ranger format
mv SRR21817808_1.fastq.gz OSA_010_S1_L002_R1_001.fastq.gz
mv SRR21817808_2.fastq.gz OSA_010_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817808

# Lane 3/4: SRR21817809
echo '  Downloading SRR21817809 (Lane 3)...'
prefetch SRR21817809 || { echo 'Failed to prefetch SRR21817809'; exit 1; }
fasterq-dump SRR21817809 --split-files --threads 8 || { echo 'Failed to dump SRR21817809'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817809_1.fastq
gzip SRR21817809_2.fastq

# Rename to Cell Ranger format
mv SRR21817809_1.fastq.gz OSA_010_S1_L003_R1_001.fastq.gz
mv SRR21817809_2.fastq.gz OSA_010_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817809

# Lane 4/4: SRR21817810
echo '  Downloading SRR21817810 (Lane 4)...'
prefetch SRR21817810 || { echo 'Failed to prefetch SRR21817810'; exit 1; }
fasterq-dump SRR21817810 --split-files --threads 8 || { echo 'Failed to dump SRR21817810'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817810_1.fastq
gzip SRR21817810_2.fastq

# Rename to Cell Ranger format
mv SRR21817810_1.fastq.gz OSA_010_S1_L004_R1_001.fastq.gz
mv SRR21817810_2.fastq.gz OSA_010_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817810

cd ../..
echo '✓ Completed OSA_010'
echo ''


#============================================================
# Sample: OSA_011 (Original: GSM6617012)
# Status: OSA, Sex: female, Age: 7
#============================================================

echo 'Processing OSA_011...'
mkdir -p fastq/OSA_011
cd fastq/OSA_011

# Lane 1/4: SRR21817799
echo '  Downloading SRR21817799 (Lane 1)...'
prefetch SRR21817799 || { echo 'Failed to prefetch SRR21817799'; exit 1; }
fasterq-dump SRR21817799 --split-files --threads 8 || { echo 'Failed to dump SRR21817799'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817799_1.fastq
gzip SRR21817799_2.fastq

# Rename to Cell Ranger format
mv SRR21817799_1.fastq.gz OSA_011_S1_L001_R1_001.fastq.gz
mv SRR21817799_2.fastq.gz OSA_011_S1_L001_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817799

# Lane 2/4: SRR21817800
echo '  Downloading SRR21817800 (Lane 2)...'
prefetch SRR21817800 || { echo 'Failed to prefetch SRR21817800'; exit 1; }
fasterq-dump SRR21817800 --split-files --threads 8 || { echo 'Failed to dump SRR21817800'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817800_1.fastq
gzip SRR21817800_2.fastq

# Rename to Cell Ranger format
mv SRR21817800_1.fastq.gz OSA_011_S1_L002_R1_001.fastq.gz
mv SRR21817800_2.fastq.gz OSA_011_S1_L002_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817800

# Lane 3/4: SRR21817801
echo '  Downloading SRR21817801 (Lane 3)...'
prefetch SRR21817801 || { echo 'Failed to prefetch SRR21817801'; exit 1; }
fasterq-dump SRR21817801 --split-files --threads 8 || { echo 'Failed to dump SRR21817801'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817801_1.fastq
gzip SRR21817801_2.fastq

# Rename to Cell Ranger format
mv SRR21817801_1.fastq.gz OSA_011_S1_L003_R1_001.fastq.gz
mv SRR21817801_2.fastq.gz OSA_011_S1_L003_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817801

# Lane 4/4: SRR21817802
echo '  Downloading SRR21817802 (Lane 4)...'
prefetch SRR21817802 || { echo 'Failed to prefetch SRR21817802'; exit 1; }
fasterq-dump SRR21817802 --split-files --threads 8 || { echo 'Failed to dump SRR21817802'; exit 1; }

# Compress FASTQ files
echo '  Compressing...'
gzip SRR21817802_1.fastq
gzip SRR21817802_2.fastq

# Rename to Cell Ranger format
mv SRR21817802_1.fastq.gz OSA_011_S1_L004_R1_001.fastq.gz
mv SRR21817802_2.fastq.gz OSA_011_S1_L004_R2_001.fastq.gz

# Clean up SRA file
rm -rf SRR21817802

cd ../..
echo '✓ Completed OSA_011'
echo ''


echo 'All samples downloaded and renamed!'
echo 'Directory structure:'
tree fastq/ -L 2
