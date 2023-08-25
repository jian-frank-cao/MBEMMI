library(parallel)
n.cores = detectCores()
library(tidyverse)
library(furrr)
plan(multisession, workers = n.cores - 1)
library(reshape2)

owncode = 5
naics_version = list(
  `v2002` = list(start = 2002, end = 2006),
  `v2007` = list(start = 2007, end = 2011),
  `v2012` = list(start = 2012, end = 2016)
)

## Read Rds files --------------------------------------------------------------
raw_deo = readRDS(
  file = "./data/raw_deo_qcew.Rds"
)


## Functions -------------------------------------------------------------------
remove_commas = function(data){
  out = gsub(",", "", data) %>%
    as.numeric
  return(out)
}


## Build deo data frames -------------------------------------------------------
data_deo = naics_version %>% 
  future_map(
    ~ {
      # locate the target years
      start = .[["start"]]
      end = .[["end"]]
      year = raw_deo %>% 
        names %>% 
        substr(., 9, 12) %>% 
        as.integer
      target =  (year >= start & year <= end) %>% 
        which
      
      # get the list of naics codes
      naics_codes = raw_deo[target] %>% 
        future_map(
          ~ {
            data = .
            out = data %>% 
              filter(Ownership_Code == owncode) %>% 
              .[c("NAICS_or_BLS_Code", "NAICS_Short_Title")] %>% 
              mutate(
                NAICS_or_BLS_Code = str_trim(NAICS_or_BLS_Code, side = "both"),
                NAICS_Short_Title = str_trim(NAICS_Short_Title, side = "both") %>% 
                  gsub("\\?", "-", .)
              ) %>% 
              unique
          }
        ) %>% 
        do.call("rbind", .) %>% 
        unique %>% 
        `rownames<-`(NULL) %>%
        group_by(NAICS_or_BLS_Code) %>% 
        mutate(
          length = str_length(NAICS_Short_Title),
          title = case_when(
            length == min(length) ~
              NAICS_Short_Title %>% as.character,
            TRUE ~ NA_character_
          )
        ) %>% 
        ungroup %>% 
        filter(!is.na(title)) %>% 
        select(c("NAICS_or_BLS_Code", "NAICS_Short_Title")) %>% 
        mutate_if(
          is.factor,
          as.character
        )
      
      # build monthly data frames
      emplvl_monthly = raw_deo[target] %>% 
        set_names(year[target]) %>% 
        future_map(
          ~ {
            data = .
            out = data %>% 
              filter(Ownership_Code == owncode) %>% 
              mutate(
                NAICS_or_BLS_Code = str_trim(NAICS_or_BLS_Code, side = "both")
              ) %>% 
              .[9:20] %>% 
              sapply(., remove_commas) %>% 
              data.frame %>% 
              .[match(naics_codes$NAICS_or_BLS_Code, data$NAICS_or_BLS_Code),] %>% 
              `colnames<-`(
                colnames(.) %>% 
                  substr(., 1, 7)
              )
          }
        ) %>% 
        do.call("cbind", .) %>% 
        cbind(naics_codes, .)

      # build qtrly data frames
      emplvl_qtrly = raw_deo[target] %>% 
        set_names(year[target]) %>% 
        future_map(
          ~ {
            data = .
            data = data %>% 
              mutate(
                NAICS_or_BLS_Code = str_trim(NAICS_or_BLS_Code, side = "both")
              ) %>% 
              filter(Ownership_Code == owncode)
            out = cbind(
              data["NAICS_or_BLS_Code"],
              data[c("Month_1_Employment", "Month_2_Employment", "Month_3_Employment")] %>% 
                sapply(., remove_commas) %>% 
                rowMeans,
              data[c("Month_4_Employment", "Month_5_Employment", "Month_6_Employment")] %>% 
                sapply(., remove_commas) %>% 
                rowMeans,
              data[c("Month_7_Employment", "Month_8_Employment", "Month_9_Employment")] %>% 
                sapply(., remove_commas) %>% 
                rowMeans,
              data[c("Month_10_Employment", "Month_11_Employment", "Month_12_Employment")] %>% 
                sapply(., remove_commas) %>% 
                rowMeans
            ) %>% 
              `colnames<-`(c("NAICS_or_BLS_Code", "qtr1", "qtr2", "qtr3", "qtr4"))
            out = out[match(naics_codes$NAICS_or_BLS_Code, out$NAICS_or_BLS_Code),] %>% 
              .[-1]
          }
        ) %>% 
        do.call("cbind", .) %>% 
        cbind(naics_codes, .)
      
      # build annual data frames
      var_list = c("Confidentiality_Flag", "Units", "Total_Wages", "Average_Monthly_Employment",
                   "Average_Annual_Wage", "Average_Weekly_Wage")
      
      data_annual = var_list %>% 
        as.list %>% 
        set_names(var_list) %>% 
        future_map(
          ~ {
            var_name = .
            out = raw_deo[target] %>% 
              set_names(year[target]) %>% 
              future_map(
                ~ {
                  data = .
                  out = data %>% 
                    mutate(
                      NAICS_or_BLS_Code = str_trim(NAICS_or_BLS_Code, side = "both")
                    ) %>% 
                    filter(Ownership_Code == owncode) %>% 
                    select(c("NAICS_or_BLS_Code", all_of(var_name)))
                  out = out[match(naics_codes$NAICS_or_BLS_Code, out$NAICS_or_BLS_Code),] %>% 
                    .[[var_name]]
                  if (var_name != "Confidentiality_Flag") {
                    out = gsub(",", "", out) %>%
                      as.numeric %>% 
                      data.frame
                  }else{
                    out = out %>% 
                      data.frame %>% 
                      mutate_if(
                        is.factor,
                        as.character
                      )
                  }
                }
              ) %>% 
              do.call("cbind", .) %>% 
              `colnames<-`(year[target])
            out = cbind(
              naics_codes,
              out
            )
          }
        )
      
      out = list(
        start = start,
        end = end,
        naics_codes = naics_codes,
        monthly = list(
          emplvl_monthly = emplvl_monthly
        ),
        qtrly = list(
          emplvl_qtrly = emplvl_qtrly
        ),
        annual = data_annual
      )
    }
  )


## Save data files -------------------------------------------------------------
saveRDS(
  data_deo,
  file = "./data/data_deo_qcew.Rds"
)
