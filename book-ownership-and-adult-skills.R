# ==============================================================================
# 청소년기 가정 내 도서 보유 수가 성인기 문해력에 미치는 영향
# PIAAC 2주기 한국 자료 분석 코드
#
# Repository : https://github.com/sangyunlee-education/book-ownership-and-adult-skills
# Data       : PIAAC Cycle 2 Korea public use file (prgkorp2.csv, OECD)
#              자료는 저장소에 포함하지 않았다. OECD에서 내려받아
#              이 스크립트와 같은 폴더에 두면 그대로 실행된다.
#
# 산출물
#   표 1. 14세 시점 도서 보유 수 범주별 배경변수의 가중 분포 및 교차분석
#   표 2. 14세 시점 도서 보유 수와 성인기 문해력의 회귀분석
#   표 3. 미관측 교란요인에 대한 민감도 분석
#
# 분석 원칙
#   - 결과변수는 문해력 유의측정값(PVLIT1-PVLIT10)을 사용한다.
#   - 최종가중치(SPFWT0)와 80개 반복가중치(SPFWT1-SPFWT80)를 적용한다.
#   - 분산 추정은 Fay의 균형반복복제법(BRR)이며 rho = 0.3을 사용한다.
#   - PV별 복합표본 추정치를 Rubin의 결합규칙으로 통합한다.
# ==============================================================================


# ==============================================================================
# 0. 실행 옵션
# ==============================================================================

DATA_PATH   <- "prgkorp2.csv"
USE_VIEWER  <- interactive()   # RStudio 표 창 자동 실행 여부
SAVE_OUTPUT <- FALSE           # TRUE이면 표를 CSV로 저장
OUTPUT_DIR  <- "output"

# 민감도 분석 대상: 기준집단(10권 이하)과 비교되는 모든 도서 보유 수 범주
SENSE_TREAT_TERMS <- c(
  "book_f11-25권",
  "book_f26-100권",
  "book_f101-200권",
  "book_f201-500권",
  "book_f500권 초과"
)
SENSE_BENCHMARK_COV <- "parent_edu_f대졸 이상"   # 기준 공변량
SENSE_KD <- c(1, 2, 3)                           # 기준 공변량 대비 배수 시나리오


# ==============================================================================
# 1. 패키지
# ==============================================================================

packages <- c("tidyverse", "survey", "sensemakr")

to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

invisible(lapply(packages, library, character.only = TRUE))

options(survey.lonely.psu = "adjust")


# ==============================================================================
# 2. 변수 설정
# ==============================================================================

pv_lit      <- paste0("PVLIT", 1:10)   # 문해력 유의측정값
main_weight <- "SPFWT0"                # 최종가중치
rep_weights <- paste0("SPFWT", 1:80)   # 반복가중치

book_var <- "J2_Q06"      # 14세 시점 가정 내 도서 보유 수
age_var  <- "AGEG5LFS"    # 연령 범주

book_levels <- c(
  "10권 이하", "11-25권", "26-100권",
  "101-200권", "201-500권", "500권 초과"
)

book_col_labels <- c(
  "10권 이하"  = "≤10",
  "11-25권"    = "11-25",
  "26-100권"   = "26-100",
  "101-200권"  = "101-200",
  "201-500권"  = "201-500",
  "500권 초과" = ">500"
)

# 통제변수(모두 범주형)
control_vars <- c(
  "gender_f",       # 성별
  "age_f",          # 연령
  "parent_edu_f",   # 부모 교육수준
  "mother_work_f",  # 14세 당시 모 경제활동
  "father_work_f",  # 14세 당시 부 경제활동
  "residence_f",    # 14세 당시 거주지역 규모
  "family14_f"      # 14세 당시 가족구조
)

scenario_levels <- paste0(SENSE_KD, "배")


# ==============================================================================
# 3. 유틸리티 함수
# ==============================================================================

as_code <- function(x) suppressWarnings(as.integer(as.character(x)))

# 문자형 특수결측 처리
recode_char_missing <- function(x) {
  if (is.character(x)) x[x %in% c(".", ".n", ".r", ".d", ".v", "", "NA")] <- NA_character_
  x
}

# 숫자형 특수결측(모름·무응답·비해당 등) 처리
recode_numeric_special_missing <- function(x) {
  x_num <- as_code(x)
  x_num[x_num %in% c(7, 8, 9, 96, 97, 98, 99, 996, 997, 998, 999)] <- NA_integer_
  x_num
}

# Fay BRR 복합표본 설계
make_design <- function(data) {
  survey::svrepdesign(
    weights          = as.formula(paste0("~", main_weight)),
    repweights       = data[, rep_weights],
    data             = data,
    type             = "Fay",
    rho              = 0.3,
    combined.weights = TRUE
  )
}

p_format <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "<.001",
    TRUE      ~ sub("^0", "", sprintf("%.3f", p))
  )
}

stars_from_p <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.10  ~ "†",
    TRUE      ~ ""
  )
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

fmt_p_with_stars <- function(p) paste0(p_format(p), stars_from_p(p))

fmt_ci <- function(low, high) paste0("[", fmt_num(low, 2), ", ", fmt_num(high, 2), "]")

show_table <- function(x, table_name, file_name = NULL) {
  cat("\n\n", table_name, "\n", sep = "")
  print(x, n = Inf, width = Inf)
  if (isTRUE(USE_VIEWER)) View(x, title = table_name)
  if (isTRUE(SAVE_OUTPUT) && !is.null(file_name)) {
    dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
    readr::write_excel_csv(x, file.path(OUTPUT_DIR, file_name))
  }
  invisible(x)
}

# 회귀식: 문해력 PV ~ 도서 보유 수 + 통제변수
make_reg_formula <- function(outcome, controls = control_vars) {
  as.formula(paste(outcome, "~ book_f +", paste(controls, collapse = " + ")))
}


# ==============================================================================
# 4. Rubin 결합규칙
# ==============================================================================

# svyglm 모형 목록(PV별)의 회귀계수 통합
pool_svyglm_models <- function(model_list, df_com, outcome_name = "문해력", n_used) {

  coef_names <- Reduce(intersect, lapply(model_list, function(x) names(coef(x))))
  if (length(coef_names) == 0) {
    stop("결합할 공통 회귀계수가 없습니다. 범주 수준 또는 모형식을 확인하세요.")
  }

  coef_mat <- do.call(rbind, lapply(model_list, function(x) coef(x)[coef_names]))
  colnames(coef_mat) <- coef_names

  vcov_list <- lapply(model_list, function(x) {
    as.matrix(vcov(x))[coef_names, coef_names, drop = FALSE]
  })

  k     <- length(model_list)
  q_bar <- colMeans(coef_mat)                 # 결합 점추정치
  u_bar <- Reduce("+", vcov_list) / k         # 표본내 분산

  if (k > 1) {
    b_mat <- stats::cov(coef_mat)             # 표본간 분산
    if (is.null(dim(b_mat))) {
      b_mat <- matrix(b_mat, nrow = 1, dimnames = list(coef_names, coef_names))
    }
  } else {
    b_mat <- matrix(0, length(q_bar), length(q_bar),
                    dimnames = list(coef_names, coef_names))
  }

  total_var_mat <- u_bar + (1 + 1 / k) * b_mat
  std_error     <- sqrt(diag(total_var_mat))

  b_diag <- diag(b_mat)
  lambda <- ((1 + 1 / k) * b_diag) / diag(total_var_mat)
  lambda <- pmin(pmax(lambda, 0), 0.999999)

  df_old    <- (k - 1) / (lambda^2)
  df_obs    <- ((df_com + 1) / (df_com + 3)) * df_com * (1 - lambda)
  df_pooled <- 1 / ((1 / df_old) + (1 / df_obs))
  df_pooled[is.na(df_pooled) | is.infinite(df_pooled)] <- df_com
  df_pooled[b_diag < 1e-12] <- df_com

  statistic <- q_bar / std_error

  tibble(
    outcome        = outcome_name,
    n              = n_used,
    term           = names(q_bar),
    estimate       = as.numeric(q_bar),
    std_error      = as.numeric(std_error),
    conf_low       = as.numeric(q_bar - qt(0.975, df = df_pooled) * std_error),
    conf_high      = as.numeric(q_bar + qt(0.975, df = df_pooled) * std_error),
    statistic_type = "t",
    statistic      = as.numeric(statistic),
    df             = as.numeric(df_pooled),
    p_value        = as.numeric(2 * pt(-abs(statistic), df = df_pooled)),
    stars          = stars_from_p(p_value)
  )
}

# 단일 추정치·표준오차 형태(민감도 분석 결과 등)의 통합
pool_scalar_df <- function(df, key_cols, estimate_col, se_col,
                           df_com, outcome_name = "문해력", n_used) {

  df_valid <- df %>%
    filter(!is.na(.data[[estimate_col]]), !is.na(.data[[se_col]]))

  if (length(key_cols) == 0) {
    grouped     <- df_valid %>% mutate(.pool_key = "all") %>% group_by(.pool_key)
    select_keys <- character(0)
  } else {
    grouped     <- df_valid %>% group_by(across(all_of(key_cols)))
    select_keys <- key_cols
  }

  pooled <- grouped %>%
    summarise(
      k        = n(),
      estimate = mean(.data[[estimate_col]]),
      u_bar    = mean((.data[[se_col]])^2),
      b_var    = if (n() > 1) stats::var(.data[[estimate_col]]) else 0,
      .groups  = "drop"
    )

  if (length(key_cols) == 0) pooled <- pooled %>% select(-.pool_key)

  pooled %>%
    mutate(
      total_var      = u_bar + (1 + 1 / k) * b_var,
      std_error      = sqrt(total_var),
      lambda         = pmin(pmax(((1 + 1 / k) * b_var) / total_var, 0), 0.999999),
      df_old         = (k - 1) / (lambda^2),
      df_obs         = ((df_com + 1) / (df_com + 3)) * df_com * (1 - lambda),
      df             = 1 / ((1 / df_old) + (1 / df_obs)),
      df             = if_else(is.na(df) | is.infinite(df) | b_var < 1e-12, df_com, df),
      conf_low       = estimate - qt(0.975, df = df) * std_error,
      conf_high      = estimate + qt(0.975, df = df) * std_error,
      statistic_type = "t",
      statistic      = estimate / std_error,
      p_value        = 2 * pt(-abs(statistic), df = df),
      stars          = stars_from_p(p_value),
      outcome        = outcome_name,
      n              = n_used
    ) %>%
    select(
      outcome, n, all_of(select_keys),
      estimate, std_error, conf_low, conf_high,
      statistic_type, statistic, df, p_value, stars
    )
}


# ==============================================================================
# 5. 자료 불러오기 및 전처리
# ==============================================================================

if (!file.exists(DATA_PATH)) {
  stop("DATA_PATH에 지정한 파일이 없습니다: ", DATA_PATH)
}

dat_raw <- readr::read_delim(
  DATA_PATH,
  delim          = NULL,   # 구분자 자동 판별
  show_col_types = FALSE,
  progress       = FALSE,
  locale         = readr::locale(encoding = "UTF-8")
)

required_raw_vars <- c(
  pv_lit, main_weight, rep_weights, book_var, age_var,
  "GENDER_R", "PAREDC2", "J2_Q04d", "J2_Q05d",
  "J2_Q07_C", "J2_Q0801", "J2_Q0802"
)

missing_raw_vars <- setdiff(required_raw_vars, names(dat_raw))
if (length(missing_raw_vars) > 0) {
  stop("원자료에서 필요한 변수를 찾지 못했습니다: ",
       paste(missing_raw_vars, collapse = ", "))
}

# 숫자형 특수결측 처리 대상(원문 설문변수)
special_numeric_vars <- c(
  book_var, "PAREDC2", "J2_Q04d", "J2_Q05d",
  "J2_Q07_C", "J2_Q0801", "J2_Q0802"
)

dat <- dat_raw %>%
  mutate(across(where(is.character), recode_char_missing)) %>%
  mutate(across(all_of(special_numeric_vars), recode_numeric_special_missing)) %>%
  mutate(
    # 14세 시점 가정 내 도서 보유 수
    book_f = case_when(
      as_code(.data[[book_var]]) == 1 ~ "10권 이하",
      as_code(.data[[book_var]]) == 2 ~ "11-25권",
      as_code(.data[[book_var]]) == 3 ~ "26-100권",
      as_code(.data[[book_var]]) == 4 ~ "101-200권",
      as_code(.data[[book_var]]) == 5 ~ "201-500권",
      as_code(.data[[book_var]]) == 6 ~ "500권 초과",
      TRUE ~ NA_character_
    ),
    book_f = factor(book_f, levels = book_levels),

    # 성별
    gender_f = case_when(
      as_code(GENDER_R) == 1 ~ "남성",
      as_code(GENDER_R) == 2 ~ "여성",
      TRUE ~ NA_character_
    ),
    gender_f = factor(gender_f, levels = c("남성", "여성")),

    # 연령
    age_f = case_when(
      as_code(.data[[age_var]]) == 1  ~ "20-24세",
      as_code(.data[[age_var]]) == 2  ~ "25-29세",
      as_code(.data[[age_var]]) == 3  ~ "30-34세",
      as_code(.data[[age_var]]) == 4  ~ "35-39세",
      as_code(.data[[age_var]]) == 5  ~ "40-44세",
      as_code(.data[[age_var]]) == 6  ~ "45-49세",
      as_code(.data[[age_var]]) == 7  ~ "50-54세",
      as_code(.data[[age_var]]) == 8  ~ "55-59세",
      as_code(.data[[age_var]]) == 9  ~ "60-64세",
      as_code(.data[[age_var]]) == 10 ~ "65세 이상",
      TRUE ~ NA_character_
    ),
    age_f = factor(
      age_f,
      levels = c("20-24세", "25-29세", "30-34세", "35-39세", "40-44세",
                 "45-49세", "50-54세", "55-59세", "60-64세", "65세 이상")
    ),

    # 부모 교육수준
    parent_edu_f = case_when(
      as_code(PAREDC2) == 1 ~ "중졸 이하",
      as_code(PAREDC2) == 2 ~ "고졸/전문대졸",
      as_code(PAREDC2) == 3 ~ "대졸 이상",
      TRUE ~ NA_character_
    ),
    parent_edu_f = factor(parent_edu_f,
                          levels = c("중졸 이하", "고졸/전문대졸", "대졸 이상")),

    # 14세 당시 모 경제활동
    mother_work_f = case_when(
      as_code(J2_Q04d) == 1 ~ "유급직",
      as_code(J2_Q04d) == 2 ~ "무직/가사 등",
      TRUE ~ NA_character_
    ),
    mother_work_f = factor(mother_work_f, levels = c("유급직", "무직/가사 등")),

    # 14세 당시 부 경제활동
    father_work_f = case_when(
      as_code(J2_Q05d) == 1 ~ "유급직",
      as_code(J2_Q05d) == 2 ~ "무직/가사 등",
      TRUE ~ NA_character_
    ),
    father_work_f = factor(father_work_f, levels = c("유급직", "무직/가사 등")),

    # 14세 당시 거주지역 규모
    residence_f = case_when(
      as_code(J2_Q07_C) == 1 ~ "대도시",
      as_code(J2_Q07_C) == 2 ~ "중소도시",
      as_code(J2_Q07_C) == 3 ~ "소도시/읍면",
      as_code(J2_Q07_C) == 4 ~ "농어촌/시골",
      TRUE ~ NA_character_
    ),
    residence_f = factor(residence_f,
                         levels = c("대도시", "중소도시", "소도시/읍면", "농어촌/시골")),

    # 14세 당시 가족구조
    family14_f = case_when(
      as_code(J2_Q0801) == 1 & as_code(J2_Q0802) == 1 ~ "양부모 동거",
      as_code(J2_Q0801) == 1 & as_code(J2_Q0802) != 1 ~ "모만 동거",
      as_code(J2_Q0801) != 1 & as_code(J2_Q0802) == 1 ~ "부만 동거",
      as_code(J2_Q0801) != 1 & as_code(J2_Q0802) != 1 ~ "생부모 비동거",
      TRUE ~ NA_character_
    ),
    family14_f = factor(family14_f,
                        levels = c("양부모 동거", "모만 동거", "부만 동거", "생부모 비동거"))
  ) %>%
  # 기준범주 설정
  mutate(
    book_f       = relevel(book_f, ref = "10권 이하"),
    parent_edu_f = relevel(parent_edu_f, ref = "중졸 이하")
  )

# 완전사례 자료
cc_data <- dat %>%
  select(all_of(c(pv_lit, "book_f", control_vars, main_weight, rep_weights))) %>%
  tidyr::drop_na() %>%
  mutate(across(where(is.factor), droplevels)) %>%
  mutate(
    book_f       = relevel(book_f, ref = "10권 이하"),
    parent_edu_f = relevel(parent_edu_f, ref = "중졸 이하")
  )

cc_n   <- nrow(cc_data)
cc_des <- make_design(cc_data)
cc_df  <- survey::degf(cc_des)

cat("\n완전사례분석 표본 크기: ", cc_n, "명\n", sep = "")
cat("복합표본 설계 자유도: ", cc_df, "\n", sep = "")


# ==============================================================================
# 6. 표 1. 도서 보유 수 범주별 배경변수의 가중 분포 및 교차분석
# ==============================================================================

background_specs <- list(
  list(var = "gender_f",      label = "성별",
       levels = c("남성", "여성")),
  list(var = "age_f",         label = "연령",
       levels = c("20-24세", "25-29세", "30-34세", "35-39세", "40-44세",
                  "45-49세", "50-54세", "55-59세", "60-64세", "65세 이상")),
  list(var = "parent_edu_f",  label = "부모 교육수준",
       levels = c("중졸 이하", "고졸/전문대졸", "대졸 이상")),
  list(var = "mother_work_f", label = "모 경제활동",
       levels = c("유급직", "무직/가사 등")),
  list(var = "father_work_f", label = "부 경제활동",
       levels = c("유급직", "무직/가사 등")),
  list(var = "residence_f",   label = "거주지",
       levels = c("대도시", "중소도시", "소도시/읍면", "농어촌/시골")),
  list(var = "family14_f",    label = "가족구조",
       levels = c("양부모 동거", "모만 동거", "부만 동거", "생부모 비동거"))
)

# Rao-Scott 수정 카이제곱 검정(통계량 옵션을 순차적으로 시도)
safe_svy_chisq_p <- function(var) {
  fml <- as.formula(paste("~ book_f +", var))
  for (method in c("F", "adjWald", "Chisq")) {
    p <- tryCatch(
      suppressWarnings(as.numeric(survey::svychisq(fml, cc_des, statistic = method)$p.value)),
      error = function(e) NA_real_
    )
    if (length(p) == 1 && is.finite(p)) return(p)
  }
  NA_real_
}

# 도서 보유 수 범주별 가중 열 백분율
weighted_col_percent <- function(var, category) {
  tmp <- cc_data %>%
    filter(!is.na(book_f), !is.na(.data[[var]])) %>%
    group_by(book_f, category_tmp = as.character(.data[[var]])) %>%
    summarise(weighted_n = sum(.data[[main_weight]], na.rm = TRUE), .groups = "drop") %>%
    group_by(book_f) %>%
    mutate(percent = 100 * weighted_n / sum(weighted_n, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(category_tmp == category)

  out <- setNames(rep(NA_real_, length(book_levels)), book_levels)
  out[as.character(tmp$book_f)] <- tmp$percent
  out
}

make_pct_row <- function(label, pct, p_text = "") {
  values <- as.list(fmt_num(pct[book_levels], 1))
  names(values) <- unname(book_col_labels[book_levels])
  bind_cols(tibble(변수 = label), as_tibble(values), tibble(`pᵃ` = p_text))
}

make_table1_section <- function(spec) {
  empty_pct <- setNames(rep(NA_real_, length(book_levels)), book_levels)

  header <- make_pct_row(spec$label, empty_pct, fmt_p_with_stars(safe_svy_chisq_p(spec$var)))
  rows   <- map_dfr(spec$levels, function(cat_label) {
    make_pct_row(cat_label, weighted_col_percent(spec$var, cat_label))
  })

  bind_rows(header, rows)
}

표1_기술통계_교차분석 <- map_dfr(background_specs, make_table1_section)

show_table(
  표1_기술통계_교차분석,
  "<표 1> 14세 시점 가정 내 도서 보유 수 범주별 주요 배경변수의 가중 분포",
  "table1.csv"
)

cat("\n표 1 주. 수치는 복합표본 최종가중치를 적용한 열 백분율(%)이다. ",
    "a Rao-Scott 수정 카이제곱 검정의 p값이다.\n", sep = "")


# ==============================================================================
# 7. 표 2. 도서 보유 수와 성인기 문해력의 회귀분석
# ==============================================================================

main_model_list <- map(pv_lit, function(y) {
  survey::svyglm(make_reg_formula(y), design = cc_des)
})

main_pooled <- pool_svyglm_models(main_model_list, df_com = cc_df, n_used = cc_n)

main_effect_results <- main_pooled %>%
  filter(str_detect(term, "^book_f"))

표2_회귀분석 <- main_effect_results %>%
  transmute(
    `도서 보유 수`  = str_remove(term, "^book_f"),
    추정치          = fmt_num(estimate, 2),
    표준오차        = fmt_num(std_error, 2),
    `t값`           = fmt_num(statistic, 2),
    `p값`           = fmt_p_with_stars(p_value),
    `95% 신뢰구간`  = fmt_ci(conf_low, conf_high)
  )

show_table(
  표2_회귀분석,
  "<표 2> 14세 시점 가정 내 도서 보유 수 범주별 성인기 문해력의 추정치",
  "table2.csv"
)

cat("\n표 2 주. 표본 크기는 ", cc_n, "명이며, 기준집단은 ‘10권 이하’이다. ", sep = "")
cat("모형은 성별, 연령, 부모 교육수준, 14세 당시 모·부의 경제활동 여부, 거주지역 규모 및 가족구조를 통제하였다. ")
cat("각 계수와 표준오차는 10개의 문해력 PV에 대해 최종가중치와 80개의 반복가중치, Fay 조정계수 0.3을 적용한 BRR로 추정한 뒤 Rubin의 결합규칙에 따라 통합하였다. ")
cat("† p < .10, * p < .05, ** p < .01, *** p < .001.\n")


# ==============================================================================
# 8. 표 3. 미관측 교란요인에 대한 민감도 분석
# ==============================================================================

# sensemakr는 lm 객체를 요구한다. 기준 공변량의 부분설명력은 최종가중치를 적용한
# lm에서 산출하고, 조정 추정치와 조정 표준오차는 복합표본 추정치를 기준으로 계산한다.
lm_model_list <- map(pv_lit, function(y) {
  stats::lm(make_reg_formula(y), data = cc_data, weights = cc_data[[main_weight]])
})

benchmark_bounds_pv <- function(i, treatment_term) {
  brr <- coef(summary(main_model_list[[i]]))
  if (!(treatment_term %in% rownames(brr))) {
    stop("벤치마크 분석 대상 계수가 모형에 없습니다: ", treatment_term)
  }

  est_brr <- brr[treatment_term, "Estimate"]
  se_brr  <- brr[treatment_term, "Std. Error"]

  sense_out <- sensemakr::sensemakr(
    model                 = lm_model_list[[i]],
    treatment             = treatment_term,
    benchmark_covariates  = SENSE_BENCHMARK_COV,
    kd                    = SENSE_KD,
    q                     = 1,
    alpha                 = 0.05
  )

  bounds <- as_tibble(sense_out$bounds)
  stopifnot(nrow(bounds) == length(SENSE_KD))

  bounds %>%
    mutate(
      treatment         = treatment_term,
      pv                = pv_lit[i],
      scenario          = scenario_levels,
      Adjusted_Estimate = sensemakr::adjusted_estimate(
        estimate = est_brr, se = se_brr, dof = cc_df,
        r2dz.x = r2dz.x, r2yz.dx = r2yz.dx
      ),
      Adjusted_SE = sensemakr::adjusted_se(
        se = se_brr, dof = cc_df,
        r2dz.x = r2dz.x, r2yz.dx = r2yz.dx
      )
    ) %>%
    select(treatment, pv, scenario, r2dz.x, r2yz.dx, Adjusted_Estimate, Adjusted_SE)
}

sensitivity_bounds_all <- map_dfr(SENSE_TREAT_TERMS, function(treatment_term) {
  map_dfr(seq_along(pv_lit), benchmark_bounds_pv, treatment_term = treatment_term)
})

sensitivity_bounds_pooled <- pool_scalar_df(
  df           = sensitivity_bounds_all,
  key_cols     = c("treatment", "scenario"),
  estimate_col = "Adjusted_Estimate",
  se_col       = "Adjusted_SE",
  df_com       = cc_df,
  n_used       = cc_n
)

# 강건성 값(RV): 통합된 계수·표준오차·자유도에 기초해 산출
sensitivity_rv_pooled <- main_effect_results %>%
  filter(term %in% SENSE_TREAT_TERMS) %>%
  mutate(
    sensitivity_stats = pmap(
      list(estimate = estimate, se = std_error, dof = df, treatment = term),
      function(estimate, se, dof, treatment) {
        sensemakr::sensitivity_stats(
          estimate = estimate, se = se, dof = dof,
          treatment = treatment, q = 1, alpha = 0.05
        )
      }
    ),
    RV        = map_dbl(sensitivity_stats, ~ as.numeric(.x$rv_q[1])),
    극단적_RV = map_dbl(sensitivity_stats, ~ as.numeric(.x$rv_qa[1])),
    `도서 보유 수` = factor(str_remove(term, "^book_f"), levels = book_levels[-1])
  ) %>%
  arrange(`도서 보유 수`) %>%
  select(treatment = term, `도서 보유 수`, RV, 극단적_RV)

표3_벤치마크 <- sensitivity_bounds_pooled %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    조정치   = paste0(fmt_num(estimate, 2), stars_from_p(p_value),
                      "(", fmt_num(std_error, 2), ")")
  ) %>%
  select(treatment, scenario, 조정치) %>%
  pivot_wider(
    names_from  = scenario,
    values_from = 조정치,
    names_glue  = "{scenario} 조정치(SE)"
  )

표3_민감도분석 <- sensitivity_rv_pooled %>%
  left_join(표3_벤치마크, by = "treatment") %>%
  mutate(
    `도서 보유 수` = as.character(`도서 보유 수`),
    RV             = fmt_num(RV, 3),
    `극단적 RV`    = fmt_num(극단적_RV, 3)
  ) %>%
  select(`도서 보유 수`, RV, `극단적 RV`, ends_with("조정치(SE)"))

show_table(
  표3_민감도분석,
  "<표 3> 미관측 교란요인에 대한 도서 보유 수 효과의 민감도 분석",
  "table3.csv"
)

cat("\n표 3 주. 기준집단은 ‘10권 이하’이다. ", sep = "")
cat("RV는 미관측 교란요인이 추정된 인과효과를 통계적으로 유의하지 않게 만들기 위해 원인변수와 결과변수를 각각 설명해야 하는 최소 부분설명력을 의미한다. ")
cat("극단적 RV는 미관측 교란요인과 결과변수의 관련성에 제한을 두지 않는 극단적 조건에서, 추정된 인과효과를 통계적으로 유의하지 않게 만들기 위해 원인변수와 가져야 하는 최소 부분설명력을 의미한다. ")
cat("1배, 2배 및 3배 시나리오는 가상의 미관측 교란요인이 ‘부모 교육수준: 대졸 이상’ 범주의 관찰된 부분설명력과 동일하거나 각각 2배 및 3배인 경우를 나타낸다. ")
cat("RV와 극단적 RV, 조정 추정치와 조정 표준오차는 10개 문해력 PV와 복합표본 설계를 반영하여 추정한 뒤 Rubin의 결합규칙에 따라 통합하였다. ")
cat("*** p < .001.\n")


# ==============================================================================
# 끝
# ==============================================================================
