# MBEMMI
MBEMMI - A Fast and Accurate Multiple Imputation Method for Missing Values in Large Multi-scale Data Sets with Linear Constraints

## Code
### ./R/step0_setup.R
**Usage**: Installs required packages if necessary.\
**Output**: None


### ./R/step1a_read_bls.R
**Usage**: Reads the raw BLS QCEW data files under "./data/BLS/".\
**Output**: ./data/raw_bls_qcew.Rds\
**Format**: List of data frames that store the raw data. Each data frame is a csv file.

### ./R/step1b_read_deo.R
**Usage**: Reads the raw DEO QCEW data files under "./data/DEO/".\
**Output**: ./data/raw_deo_qcew.Rds\
**Format**: List of data frames that store the raw data. Each data frame is a csv file.

### ./R/step2a_preprocess_bls.R
**Usage**: Reads "./data/raw_bls_qcew.Rds" and reformats the data into time series.\
**Output**: ./data/data_bls_qcew.Rds\
**Format**: List of data frames that store the time series. In each data frame, rows are NAICS codes, columns are time. The data frames are grouped by NAICS versions (v2002, v2007, v2012) and frequency (monthly, quarterly, annual).

### ./R/step2b_preprocess_deo.R
**Usage**: Reads "./data/raw_deo_qcew.Rds" and reformats the data into time series.\
**Output**: ./data/data_deo_qcew.Rds\
**Format**: List of data frames that store the time series. In each data frame, rows are NAICS codes, columns are time. The data frames are grouped by NAICS versions (v2002, v2007, v2012) and frequency (monthly, quarterly, annual).

### ./R/step3_merge.R
**Usage**: Reads "./data/data_bls_qcew.Rds" and "./data/data_deo_qcew.Rds" and constructs new employment time series using unsuppressed employment data from DEO and suppression positions from BLS. It also separates the mixed 2-digit industries (e.g. seperates industry "31-33" into 31, 32, and 33).\
**Output**: ./data/qcew_2digit_fixed.Rds\
**Format**: List of data frames that store the time series. In each data frame, row is NAICS code, column is time. It contains BLS suppressed DEO employment data, unsuppressed DEO employment data, and unsuppressed BLS establishment data.

### ./R/step4_suppression.R
**Usage**: Reads "./data/qcew_2digit_fixed.Rds" and applies random primary suppression and recursive secondary suppression. Outputs 50 randomly suppressed data sets. \
**Output**: ./data/qcew_rnd_suppression.Rds\
**Format**: List of 50 randomly suppressed data sets. Each data set has the same format as in "./data/qcew_2digit_fixed.Rds".

### ./R/step5_imputatoin.R (need to be fixed)
**Usage**: Reads "./data/qcew_rnd_suppression.Rds" and applies PSI-MBEMMI, BMMI, and EMB method to obtain imputed data sets. \
**Output**: ./data/QCEW_imputations.Rds\
**Format**: List of N imputed data frames. In each data frame, rows are NAICS codes, columns are time.®
