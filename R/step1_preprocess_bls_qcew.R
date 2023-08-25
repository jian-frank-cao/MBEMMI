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
raw_bls = readRDS(
  file = "./data/raw_bls_qcew.Rds"
)


## Build bls data frames -------------------------------------------------------
data_bls = naics_version %>% 
  future_map(
    ~ {
      # locate the target years
      start = .[["start"]]
      end = .[["end"]]
      year = raw_bls$raw_bls_qtrly %>% 
        names %>% 
        substr(., 1, 4) %>% 
        as.integer
      target =  (year >= start & year <= end) %>% 
        which
      
      # get the list of naics codes
      naics_codes = raw_bls$raw_bls_qtrly[target] %>% 
        future_map(
          ~ {
            data = .
            out = data %>% 
              filter(own_code == owncode) %>% 
              .[c("industry_code", "industry_title")] %>% 
              unique
          }
        ) %>% 
        do.call("rbind", .) %>% 
        unique %>% 
        `rownames<-`(NULL) %>%
        group_by(industry_code) %>% 
        mutate(
          length = str_length(industry_title),
          title = case_when(
            length == min(length) ~
              industry_title %>% as.character,
            TRUE ~ NA_character_
          )
        ) %>% 
        ungroup %>% 
        filter(!is.na(title)) %>% 
        select(c("industry_code", "industry_title")) %>% 
        mutate_if(
          is.factor,
          as.character
        )
      
      # build quarterly and monthly data frames
      var_list_detail = c("disclosure_code", "qtrly_estabs_count",
                          "emplvl_qtrly", "emplvl_monthly",
                          "total_qtrly_wages", "avg_wkly_wage")
      
      data_detail = var_list_detail %>% 
        as.list %>% 
        set_names(var_list_detail) %>% 
        future_map(
          ~ {
            var_name = .
            out = raw_bls$raw_bls_qtrly[target] %>% 
              set_names(year[target]) %>% 
              future_map(
                ~ {
                  data = .
                  if (var_name == "emplvl_qtrly") {
                    out = data %>% 
                      filter(own_code == owncode) %>% 
                      melt(
                        .,
                        id.vars = c("industry_code", "qtr"),
                        measure.vars = c("month1_emplvl", "month2_emplvl",
                                         "month3_emplvl")
                      ) %>% 
                      dcast(
                        .,
                        industry_code ~ qtr,
                        mean
                      )%>% 
                      `colnames<-`(c("industry_code", "qtr1", "qtr2",
                                     "qtr3", "qtr4"))
                  }else if(var_name == "emplvl_monthly"){
                    out = data %>% 
                      filter(own_code == owncode) %>% 
                      melt(
                        .,
                        id.vars = c("industry_code", "qtr"),
                        measure.vars = c("month1_emplvl", "month2_emplvl",
                                         "month3_emplvl")
                      ) %>% 
                      dcast(
                        .,
                        industry_code ~ qtr + variable
                      )
                    colnames(out)[2:13] = paste0(
                      "qtr",
                      colnames(out)[2:13] %>%
                        substr(., 1, 8)
                    )
                  }else{
                    out = data %>% 
                      filter(own_code == owncode) %>% 
                      dcast(
                        .,
                        industry_code ~ qtr,
                        value.var = var_name
                      ) %>% 
                      `colnames<-`(c("industry_code", "qtr1", "qtr2",
                                     "qtr3", "qtr4"))
                  }
                  
                  out = out[match(naics_codes$industry_code, out$industry_code),
                            2:ncol(out)]
                }
              ) %>% 
              do.call("cbind", .)
            out = cbind(
              naics_codes,
              out
            ) %>% 
              mutate_if(
                is.factor,
                as.character
              )
          }
        )
      
      # build annual data frames
      var_list_annual = c("disclosure_code", "annual_avg_estabs_count",
                          "annual_avg_emplvl", "total_annual_wages",
                          "annual_avg_wkly_wage")
      
      data_annual = var_list_annual %>% 
        as.list %>% 
        set_names(var_list_annual) %>% 
        future_map(
          ~ {
            var_name = .
            out = raw_bls$raw_bls_annual[target] %>% 
              future_map(
                ~ {
                  data = .
                  out = data %>% 
                    filter(own_code == owncode) %>% 
                    select(c("industry_code", all_of(var_name)))
                  out = out[match(naics_codes$industry_code,
                                  out$industry_code),] %>% 
                    select(all_of(var_name))
                }
              ) %>% 
              do.call("cbind", .) %>% 
              `colnames<-`(year[target])
            out = cbind(
              naics_codes,
              out
            ) %>% 
              mutate_if(
                is.factor,
                as.character
              )
          }
        )
      
      out = list(
        start = start,
        end = end,
        naics_codes = naics_codes,
        monthly = data_detail["emplvl_monthly"],
        qtrly = data_detail[c("disclosure_code", "qtrly_estabs_count",
                              "emplvl_qtrly", "total_qtrly_wages",
                              "avg_wkly_wage")],
        annual = data_annual
      )
    }
  )


## Save data files -------------------------------------------------------------
saveRDS(
  data_bls,
  file = "./data/data_bls_qcew.Rds"
)
