################################################################################
# Replication Code for:
# "청소년기 가정 내 도서보유수가 성인기 정보처리스킬에 미치는 영향"
#
# Data: PIAAC Cycle 2 (Korea)
################################################################################

# ==============================================================================
# 0. 패키지 및 초기 설정
# ==============================================================================
library(tidyverse)
library(survey)
library(mitools)
library(readr)
library(stringr)
library(sensemakr)
library(emmeans)
library(ggplot2)
library(showtext)

# 복합표본 PSU 조정 옵션
options(survey.lonely.psu = "adjust")

# 그래프 폰트 및 해상도 설정
font_add_google("Nanum Gothic", "nanumgothic")
showtext_auto()
showtext_opts(dpi = 600)

# emmeans 조합 제한 해제
emm_options(rg.limit = 500000)

# ==============================================================================
# 1. 데이터 불러오기 및 전처리
# ==============================================================================
# 데이터 로드 (경로는 작업 디렉토리에 맞게 수정)
dat_raw <- read_delim(
  "prgkorp2.csv",
  delim = ";",
  show_col_types = FALSE,
  locale = locale(encoding = "UTF-8")
)

# 주요 변수명 설정
pv_lit <- paste0("PVLIT", 1:10)   # 언어능력 PVs
pv_num <- paste0("PVNUM", 1:10)   # 수리력 PVs
pv_aps <- paste0("PVAPS", 1:10)   # 적응적 문제해결력 PVs

main_weight <- "SPFWT0"
rep_weights <- paste0("SPFWT", 1:80)
book_var <- "J2_Q06"
age_var <- "AGEG5LFS"

base_controls <- c("GENDER_R", age_var, "PAREDC2")
extended_controls <- c("GENDER_R", age_var, "PAREDC2", "J2_Q04d", "J2_Q05d", "J2_Q07_C", "family14")

# 존재하는 변수만 필터링하는 함수
keep_existing <- function(vars, data) {
  vars[vars %in% names(data)]
}

# 특수결측치 처리 함수 (.n, .d 등)
recode_missing <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    x[x %in% c(7, 8, 9, 96, 97, 98, 99, 996, 997, 998, 999)] <- NA
  }
  if (is.character(x)) {
    x[x %in% c(".n", ".r", ".d", ".v", "", "NA")] <- NA
  }
  return(x)
}

# 결측치 정제 적용
dat <- dat_raw %>%
  mutate(
    across(
      all_of(keep_existing(unique(c(book_var, extended_controls)), dat_raw)),
      recode_missing
    )
  )

# 파생변수 생성 (도서보유수 및 14 세 당시 가족구조)
dat <- dat %>%
  mutate(
    J2_Q06_f = case_when(
      .data[[book_var]] == 1 ~ "10 권 이하",
      .data[[book_var]] == 2 ~ "11-25 권",
      .data[[book_var]] == 3 ~ "26-100 권",
      .data[[book_var]] == 4 ~ "101-200 권",
      .data[[book_var]] == 5 ~ "201-500 권",
      .data[[book_var]] == 6 ~ "500 권 초과",
      TRUE ~ NA_character_
    ),
    J2_Q06_f = factor(
      J2_Q06_f,
      levels = c("10 권 이하", "11-25 권", "26-100 권", "101-200 권", "201-500 권", "500 권 초과")
    ),
    family14 = case_when(
      J2_Q0801 == 1 & J2_Q0802 == 1 ~ "양부모 동거",
      J2_Q0801 == 1 & J2_Q0802 != 1 ~ "모만 동거",
      J2_Q0801 != 1 & J2_Q0802 == 1 ~ "부만 동거",
      J2_Q0801 != 1 & J2_Q0802 != 1 ~ "생부모 비동거",
      TRUE ~ NA_character_
    ),
    family14 = factor(
      family14,
      levels = c("양부모 동거", "모만 동거", "부만 동거", "생부모 비동거")
    )
  )

# 범주형 변수 팩터 처리
categorical_vars <- keep_existing(extended_controls, dat)
dat <- dat %>%
  mutate(across(all_of(categorical_vars), as.factor)) %>%
  mutate(across(where(is.factor), droplevels))

# ==============================================================================
# 2. 복합표본 설계 (Design) 함수 정의
# ==============================================================================
make_design <- function(data) {
  svrepdesign(
    weights = as.formula(paste0("~", main_weight)),
    repweights = data[, rep_weights],
    data = data,
    type = "Fay",
    rho = 0.3, # Fay 조정계수
    combined.weights = TRUE
  )
}

# ==============================================================================
# 3. 기초 분석: 원인변수별 외생적 배경변수의 기술통계 및 교차분석 결과
# ==============================================================================
run_homogeneity_test <- function(data) {
  # 분석 대상 유효 결측치 제거
  dat_clean <- data %>%
    filter(complete.cases(select(., all_of(c("J2_Q06_f", extended_controls))))) %>%
    mutate(across(all_of(c("J2_Q06_f", extended_controls)), droplevels))
  
  des_clean <- make_design(dat_clean)
  
  target_vars <- extended_controls
  results_list <- list()
  
  for (v in target_vars) {
    fml <- as.formula(paste("~", v, "+ J2_Q06_f"))
    svy_tab <- svytable(fml, design = des_clean)
    svy_prop <- round(prop.table(svy_tab, margin = 2) * 100, 2)
    chisq_res <- svychisq(fml, design = des_clean, statistic = "F")
    
    merged <- as.data.frame(svy_tab) %>%
      inner_join(
        as.data.frame(svy_prop),
        by = c(v, "J2_Q06_f"),
        suffix = c("_Freq", "_Pct")
      ) %>%
      mutate(
        변수명 = v,
        범주 = .[, v],
        출력 = paste0(round(Freq_Freq, 0), " (", round(Freq_Pct, 1), "%)"),
        p_value = chisq_res$p.value
      ) %>%
      select(변수명, 범주, J2_Q06_f, 출력, p_value)
    
    results_list[[v]] <- merged
  }
  
  bind_rows(results_list) %>%
    pivot_wider(
      id_cols = c(변수명, 범주, p_value),
      names_from = J2_Q06_f,
      values_from = 출력
    )
}

table1_homogeneity <- run_homogeneity_test(dat)
print(table1_homogeneity, n = Inf)

# ==============================================================================
# 4. [연구문제 1] 다중회귀분석
# ==============================================================================
fit_pv_regression <- function(data, pv_vars, controls, outcome_name) {
  needed <- c(pv_vars, "J2_Q06_f", controls, main_weight, rep_weights)
  data_sub <- data %>% select(all_of(needed)) %>% drop_na()
  
  des <- make_design(data_sub)
  models <- vector("list", length(pv_vars))
  
  for (i in seq_along(pv_vars)) {
    fml <- as.formula(paste(pv_vars[i], "~ J2_Q06_f +", paste(controls, collapse = " + ")))
    models[[i]] <- svyglm(fml, design = des)
  }
  
  # Rubin 의 결합규칙 적용
  combined <- MIcombine(models)
  
  out <- as.data.frame(summary(combined)) %>%
    mutate(
      term = rownames(.),
      outcome = outcome_name,
      n = nrow(data_sub),
      estimate = results,
      # t-값 및 p-값 계산 (PIAAC 복합표본 자유도 80 기준)
      t_val = estimate / se,
      p_val = 2 * pt(-abs(t_val), df = 80),
      # p-값에 따른 유의성 별표 부여
      stars = case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        TRUE          ~ ""
      )
    ) %>%
    filter(str_detect(term, "J2_Q06_f")) %>%
    select(outcome, n, term, estimate, se, p_val, stars)
  
  return(out)
}

# 3 가지 역량 영역 결과 병합
rq1_results <- bind_rows(
  fit_pv_regression(dat, pv_lit, extended_controls, "언어능력"),
  fit_pv_regression(dat, pv_num, extended_controls, "수리력"),
  fit_pv_regression(dat, pv_aps, extended_controls, "적응적 문제해결력")
)

# 논문 표 (Table) 양식
rq1_results_formatted <- rq1_results %>%
  mutate(
    estimate_fmt = sprintf("%.2f", estimate),
    se_fmt = sprintf("%.2f", se),
    `B(SE)` = paste0(estimate_fmt, stars, " (", se_fmt, ")")
  ) %>%
  select(역량영역 = outcome, 표본크기 = n, 통제변수 = term, `B(SE)`, p_val)

# 결과 출력
print(as_tibble(rq1_results_formatted), n = Inf)

# ==============================================================================
# 5-1. [사전 분석] 민감도 분석을 위한 기준-공변량 진단
# ==============================================================================
message("\n================ [사전 진단] 통제변수별 교란 위험도 (t-값) 추출 ================")

check_data <- dat %>%
  select(all_of(c("PVLIT1", "J2_Q06_f", extended_controls, main_weight))) %>%
  drop_na() %>%
  mutate(
    D_bin = ifelse(J2_Q06_f == "201-500 권", 1, 0)
  )

fit_Y <- lm(
  PVLIT1 ~ J2_Q06_f + GENDER_R + AGEG5LFS + PAREDC2 +
    J2_Q04d + J2_Q05d + J2_Q07_C + family14,
  data = check_data,
  weights = check_data[[main_weight]]
)

coef_Y <- as.data.frame(coef(summary(fit_Y))) %>%
  mutate(Variable = rownames(.), Abs_t_Y = abs(`t value`)) %>%
  select(Variable, Abs_t_Y)

fit_D <- lm(
  D_bin ~ GENDER_R + AGEG5LFS + PAREDC2 +
    J2_Q04d + J2_Q05d + J2_Q07_C + family14,
  data = check_data,
  weights = check_data[[main_weight]]
)

coef_D <- as.data.frame(coef(summary(fit_D))) %>%
  mutate(Variable = rownames(.), Abs_t_D = abs(`t value`)) %>%
  select(Variable, Abs_t_D)

confounder_ranking <- coef_Y %>%
  inner_join(coef_D, by = "Variable") %>%
  filter(!str_detect(Variable, "Intercept|J2_Q06_f")) %>%
  mutate(
    Confounding_Risk = Abs_t_Y * Abs_t_D
  ) %>%
  arrange(desc(Confounding_Risk)) %>%
  select(
    통제변수 = Variable,
    `종속변수_t 값 (Y)` = Abs_t_Y,
    `처치변수_t 값 (D)` = Abs_t_D,
    총위험도 = Confounding_Risk
  )

confounder_ranking_formatted <- confounder_ranking %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

print(as_tibble(confounder_ranking_formatted), n = 10)

# ==============================================================================
# 5. [연구문제 2] 민감도 분석 (Cinelli & Hazlett) - 복합표본 SE 반영
# ==============================================================================
run_sensemakr_brr <- function(data, pv_vars, controls, outcome_name, treat_term, benchmark_cov) {
  needed <- c(pv_vars, "J2_Q06_f", controls, main_weight, rep_weights)
  data_sub <- data %>% select(all_of(needed)) %>% drop_na()
  
  des <- make_design(data_sub)
  dof <- degf(des)
  
  pv_bound_list <- list()
  orig_list <- list()
  
  for (i in seq_along(pv_vars)) {
    fml_svy <- as.formula(paste(pv_vars[i], "~ J2_Q06_f +", paste(controls, collapse = " + ")))
    
    fit_svy <- svyglm(fml_svy, design = des)
    svy_coef <- coef(summary(fit_svy))
    est_brr <- svy_coef[treat_term, "Estimate"]
    se_brr <- svy_coef[treat_term, "Std. Error"]
    
    fit_lm <- lm(fml_svy, data = data_sub, weights = data_sub[[main_weight]])
    sense_out <- sensemakr(
      model = fit_lm,
      treatment = treat_term,
      benchmark_covariates = benchmark_cov,
      kd = c(1, 2, 3)
    )
    bnd <- sense_out$bounds
    
    bnd$adj_est_brr <- sensemakr::adjusted_estimate(
      estimate = est_brr,
      se = se_brr,
      dof = dof,
      r2dz.x = bnd$r2dz.x,
      r2yz.dx = bnd$r2yz.dx
    )
    bnd$adj_se_brr <- sensemakr::adjusted_se(
      se = se_brr,
      dof = dof,
      r2dz.x = bnd$r2dz.x,
      r2yz.dx = bnd$r2yz.dx
    )
    bnd$pv <- i
    
    pv_bound_list[[i]] <- bnd
    orig_list[[i]] <- tibble(pv = i, orig_est = est_brr, orig_se = se_brr)
  }
  
  adj_results <- bind_rows(pv_bound_list) %>%
    group_by(bound_label) %>%
    summarise(
      adj_est_mean = mean(adj_est_brr),
      u_bar = mean(adj_se_brr^2),
      b_var = var(adj_est_brr),
      adj_se_total = sqrt(u_bar + (1 + 1/n()) * b_var),
      .groups = "drop"
    )
  
  orig_results <- bind_rows(orig_list) %>%
    summarise(
      orig_est_mean = mean(orig_est),
      u_bar = mean(orig_se^2),  # ← 수정: se^2 → orig_se^2
      b_var = var(orig_est),
      orig_se_total = sqrt(u_bar + (1 + 1/n()) * b_var)
    )
  
  adj_results %>%
    mutate(
      outcome = outcome_name,
      Original_Estimate = orig_results$orig_est_mean,
      Original_SE = orig_results$orig_se_total,
      Adjusted_Estimate = adj_est_mean,
      Adjusted_SE = adj_se_total
    ) %>%
    select(outcome, bound_label, Original_Estimate, Original_SE, Adjusted_Estimate, Adjusted_SE)
}

# 실행
rq2_bounds <- bind_rows(
  run_sensemakr_brr(dat, pv_lit, extended_controls, "언어능력", "J2_Q06_f201-500 권", "PAREDC23"),
  run_sensemakr_brr(dat, pv_num, extended_controls, "수리력", "J2_Q06_f201-500 권", "PAREDC23"),
  run_sensemakr_brr(dat, pv_aps, extended_controls, "적응적 문제해결력", "J2_Q06_f201-500 권", "PAREDC23")
)

# 논문 표 (Table) 양식
rq2_bounds_formatted <- rq2_bounds %>%
  mutate(
    t_val = Adjusted_Estimate / Adjusted_SE,
    p_val = 2 * pt(-abs(t_val), df = 80),
    stars = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      p_val < 0.10  ~ "†",
      TRUE          ~ ""
    ),
    Adj_Est_fmt = sprintf("%.2f", Adjusted_Estimate),
    Adj_SE_fmt = sprintf("%.2f", Adjusted_SE),
    `Badj(SEadj)` = paste0(Adj_Est_fmt, stars, " (", Adj_SE_fmt, ")")
  ) %>%
  select(
    역량영역 = outcome,
    시나리오 = bound_label,
    `Badj(SEadj)`,
    p_val
  )

print(as_tibble(rq2_bounds_formatted), n = Inf)


# ==============================================================================
# 6. [연구문제 3] 상호작용 및 사후 대비 분석
# ==============================================================================
controls_no_pared <- setdiff(extended_controls, "PAREDC2")

# 6-1. 상호작용항 통계적 유의성 검정
message("\n=== [연구문제 3-1] 상호작용항 검정 결과 ===")
fit_interaction_list <- list()
for (v in c(pv_lit, pv_num, pv_aps)) {
  fml <- as.formula(paste(v, "~ J2_Q06_f * PAREDC2 +", paste(controls_no_pared, collapse = " + ")))
  fit_interaction_list[[v]] <- svyglm(fml, design = make_design(dat))
}

# MIcombine 으로 결합하여 상호작용항 통계량 추출
interaction_summary <- summary(MIcombine(fit_interaction_list))
interaction_significant <- as_tibble(interaction_summary, rownames = "term") %>%
  filter(str_detect(term, ":"))

print(interaction_significant)

# 6-2. 추정 평균 (Estimated Marginal Means) 시각화 함수
plot_emmeans <- function(pv_vars, target_outcome, file_suffix) {
  needed <- c(pv_vars, "J2_Q06_f", "PAREDC2", controls_no_pared, main_weight, rep_weights)
  data_sub <- dat %>%
    select(all_of(needed)) %>%
    drop_na() %>%
    mutate(across(all_of(c("J2_Q06_f", "PAREDC2", controls_no_pared)), droplevels))
  
  des <- make_design(data_sub)
  emm_list <- list()
  
  for (i in seq_along(pv_vars)) {
    fml <- as.formula(paste(pv_vars[i], "~ J2_Q06_f * PAREDC2 +", paste(controls_no_pared, collapse = " + ")))
    fit <- svyglm(fml, design = des)
    fit$data <- des$variables
    emm_res <- as.data.frame(emmeans(fit, ~ J2_Q06_f | PAREDC2, data = des$variables))
    emm_res$pv <- i
    emm_list[[i]] <- emm_res
  }
  
  plot_data <- bind_rows(emm_list) %>%
    group_by(PAREDC2, J2_Q06_f) %>%
    summarise(mean_score = mean(emmean), .groups = "drop") %>%
    mutate(
      parent_edu = factor(
        case_when(
          PAREDC2 == 1 ~ "부모 중졸 이하",
          PAREDC2 == 2 ~ "부모 고졸/전문대졸",
          PAREDC2 == 3 ~ "부모 대졸 이상"
        ),
        levels = c("부모 중졸 이하", "부모 고졸/전문대졸", "부모 대졸 이상")
      ),
      J2_Q06_f = factor(
        J2_Q06_f,
        levels = c("10 권 이하", "11-25 권", "26-100 권", "101-200 권", "201-500 권", "500 권 초과")
      )
    )
  
  p <- ggplot(plot_data, aes(x = J2_Q06_f, y = mean_score, group = parent_edu, color = parent_edu, fill = parent_edu)) +
    geom_line(linewidth = 1.5, position = position_dodge(width = 0.15)) +
    geom_point(size = 4, shape = 21, color = "white", stroke = 1.2, position = position_dodge(width = 0.15)) +
    scale_color_manual(values = c("#D35400", "#2980B9", "#27AE60")) +
    scale_fill_manual(values = c("#D35400", "#2980B9", "#27AE60")) +
    theme_minimal(base_family = "nanumgothic") +
    labs(
      x = "14 세 시점 가정 내 도서보유수",
      y = paste0(target_outcome, " 추정 평균 점수"),
      color = "부모 교육수준",
      fill = "부모 교육수준"
    ) +
    theme(legend.position = "bottom")
  
  ggsave(
    filename = paste0("Figure_", file_suffix, ".png"),
    plot = p,
    width = 10,
    height = 7,
    dpi = 600
  )
  return(p)
}

p_lit <- plot_emmeans(pv_lit, "언어능력", "literacy")
p_num <- plot_emmeans(pv_num, "수리력", "numeracy")
p_aps <- plot_emmeans(pv_aps, "적응적 문제해결력", "aps")

# 6-3. 사후 대비 분석
run_contrast_analysis <- function(pv_vars, outcome_name) {
  needed <- c(pv_vars, "J2_Q06_f", "PAREDC2", controls_no_pared, main_weight, rep_weights)
  data_sub <- dat %>%
    select(all_of(needed)) %>%
    drop_na() %>%
    mutate(across(all_of(c("J2_Q06_f", "PAREDC2", controls_no_pared)), droplevels))
  
  des <- make_design(data_sub)
  
  contrast_types <- c(
    "vs_대졸이상_도서10권이하",
    "vs_고졸전문대졸_도서10권이하",
    "vs_대졸이상_도서11-25권",
    "vs_고졸전문대졸_도서11-25권"
  )
  
  est_mat <- matrix(
    0,
    nrow = length(pv_vars),
    ncol = 4,
    dimnames = list(NULL, contrast_types)
  )
  se_mat <- matrix(
    0,
    nrow = length(pv_vars),
    ncol = 4,
    dimnames = list(NULL, contrast_types)
  )
  
  for (i in seq_along(pv_vars)) {
    fml <- as.formula(paste(pv_vars[i], "~ J2_Q06_f * PAREDC2 +", paste(controls_no_pared, collapse = " + ")))
    fit <- svyglm(fml, design = des)
    fit$data <- des$variables
    
    emm <- emmeans(fit, ~ J2_Q06_f * PAREDC2, data = des$variables)
    emm_df <- as.data.frame(emm)
    
    idx_target <- which(emm_df$J2_Q06_f == "201-500 권" & emm_df$PAREDC2 == "1")
    idx_comp_A <- which(emm_df$J2_Q06_f == "10 권 이하" & emm_df$PAREDC2 == "3")
    idx_comp_B <- which(emm_df$J2_Q06_f == "10 권 이하" & emm_df$PAREDC2 == "2")
    idx_comp_C <- which(emm_df$J2_Q06_f == "11-25 권" & emm_df$PAREDC2 == "3")
    idx_comp_D <- which(emm_df$J2_Q06_f == "11-25 권" & emm_df$PAREDC2 == "2")
    
    vec_A <- rep(0, nrow(emm_df)); vec_A[idx_target] <- 1; vec_A[idx_comp_A] <- -1
    vec_B <- rep(0, nrow(emm_df)); vec_B[idx_target] <- 1; vec_B[idx_comp_B] <- -1
    vec_C <- rep(0, nrow(emm_df)); vec_C[idx_target] <- 1; vec_C[idx_comp_C] <- -1
    vec_D <- rep(0, nrow(emm_df)); vec_D[idx_target] <- 1; vec_D[idx_comp_D] <- -1
    
    contrast_res <- as.data.frame(
      contrast(emm, method = list("vs_A" = vec_A, "vs_B" = vec_B, "vs_C" = vec_C, "vs_D" = vec_D))
    )
    
    est_mat[i, ] <- contrast_res$estimate
    se_mat[i, ] <- contrast_res$SE
  }
  
  results <- list()
  for (type in contrast_types) {
    m_est <- mean(est_mat[, type])
    w_var <- mean(se_mat[, type]^2)
    b_var <- var(est_mat[, type])
    t_se <- sqrt(w_var + (1 + 1/length(pv_vars)) * b_var)
    t_val <- m_est / t_se
    p_val <- 2 * (1 - pnorm(abs(t_val)))
    
    results[[type]] <- tibble(
      대비유형 = type,
      점수차이_B = m_est,
      표준오차_SE = t_se,
      t값 = t_val,
      p값 = p_val
    )
  }
  
  bind_rows(results) %>%
    mutate(역량영역 = outcome_name) %>%
    select(역량영역, everything())
}

rq3_contrasts <- bind_rows(
  run_contrast_analysis(pv_lit, "언어능력"),
  run_contrast_analysis(pv_num, "수리력"),
  run_contrast_analysis(pv_aps, "적응적 문제해결력")
)

# 논문 표 (Table) 양식
rq3_contrasts_formatted <- rq3_contrasts %>%
  mutate(
    stars = case_when(
      p값 < 0.001 ~ "***",
      p값 < 0.01  ~ "**",
      p값 < 0.05  ~ "*",
      p값 < 0.10  ~ "†",
      TRUE        ~ ""
    ),
    B_fmt = sprintf("%.2f", 점수차이_B),
    SE_fmt = sprintf("%.2f", 표준오차_SE),
    `B(SE)` = paste0(B_fmt, stars, " (", SE_fmt, ")")
  ) %>%
  select(
    역량영역,
    대비유형,
    `B(SE)`,
    p값
  )

print(as_tibble(rq3_contrasts_formatted), n = Inf)
