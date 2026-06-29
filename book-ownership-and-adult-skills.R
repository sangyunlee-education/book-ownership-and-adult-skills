################################################################################
# 청소년기 가정 내 도서보유수와 성인기 문해력 분석 R 코드
#
# 분석 순서:
#   표 1. 14세 시점 도서 보유 수별 외생적 배경변수 기술통계 및 교차분석
#   표 2. 14세 시점 도서 보유 수와 성인기 문해력의 회귀분석
#   표 3. 미관측 교란요인에 대한 민감도 분석
#   표 4. 도서 보유 수와 부모 교육수준의 상호작용항
#   표 5. 도서풍부 취약계층과 도서빈약 상위배경 집단의 사후 대비
#
# 분석 원칙:
#   - 종속변수는 문해력 plausible values(PVLIT1-PVLIT10)를 사용한다.
#   - 최종가중치(SPFWT0)와 80개 반복가중치(SPFWT1-SPFWT80)를 적용한다.
#   - 반복가중치 방식은 Fay의 BRR이며 rho = 0.3을 사용한다.
#   - 각 PV별 복합표본 회귀 결과를 추정한 뒤 Rubin 결합규칙으로 통합한다.
################################################################################


# ==============================================================================
# 0. 실행 옵션
# ==============================================================================

DATA_PATH <- "prgkorp2.csv"

# RStudio에서 표 창을 자동으로 띄울지 여부
USE_VIEWER <- interactive()

# 민감도 분석 대상 효과
SENSE_TREAT_TERM <- "book_f201-500권"
SENSE_BENCHMARK_COV <- "parent_edu_f대졸 이상"
SENSE_KD <- c(1, 2, 3)

# 사후대비 유의확률 보정 방식
POSTHOC_ADJUST_METHOD <- "holm"


# ==============================================================================
# 1. 패키지 설치 및 로드
# ==============================================================================

packages <- c(
  "tidyverse",
  "survey",
  "sensemakr",
  "emmeans"
)

new_packages <- packages[!(packages %in% rownames(installed.packages()))]
if (length(new_packages) > 0) {
  install.packages(new_packages, dependencies = TRUE)
}

invisible(lapply(packages, library, character.only = TRUE))

options(survey.lonely.psu = "adjust")
emmeans::emm_options(rg.limit = 500000)


# ==============================================================================
# 2. 주요 변수 설정
# ==============================================================================

# 문해력 plausible values
pv_lit <- paste0("PVLIT", 1:10)

# 가중치
main_weight <- "SPFWT0"
rep_weights <- paste0("SPFWT", 1:80)

# 원자료 변수명
book_var <- "J2_Q06"
age_var <- "AGEG5LFS"

# 도서 보유 수 범주
book_levels <- c(
  "10권 이하",
  "11-25권",
  "26-100권",
  "101-200권",
  "201-500권",
  "500권 초과"
)

book_col_labels <- c(
  "10권 이하" = "≤10",
  "11-25권" = "11-25",
  "26-100권" = "26-100",
  "101-200권" = "101-200",
  "201-500권" = "201-500",
  "500권 초과" = ">500"
)

# 회귀분석 통제변수: 모두 범주형
control_vars <- c(
  "gender_f",
  "age_f",
  "parent_edu_f",
  "mother_work_f",
  "father_work_f",
  "residence_f",
  "family14_f"
)

# 상호작용 모형에서는 부모 교육수준을 상호작용항에 포함하므로 통제변수에서 제외
control_vars_no_parent_edu <- setdiff(control_vars, "parent_edu_f")


# ==============================================================================
# 3. 유틸리티 함수
# ==============================================================================

as_code <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

recode_char_missing <- function(x) {
  if (is.character(x)) {
    x[x %in% c(".", ".n", ".r", ".d", ".v", "", "NA")] <- NA_character_
  }
  x
}

recode_numeric_special_missing <- function(x) {
  x_num <- as_code(x)
  x_num[x_num %in% c(7, 8, 9, 96, 97, 98, 99, 996, 997, 998, 999)] <- NA_integer_
  x_num
}

make_design <- function(data) {
  survey::svrepdesign(
    weights = as.formula(paste0("~", main_weight)),
    repweights = data[, rep_weights],
    data = data,
    type = "Fay",
    rho = 0.3,
    combined.weights = TRUE
  )
}

p_format <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "< .001",
    TRUE ~ sub("^0", "", sprintf("%.3f", p))
  )
}

stars_from_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.10  ~ "†",
    TRUE       ~ ""
  )
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

fmt_p_with_stars <- function(p) {
  paste0(p_format(p), stars_from_p(p))
}

fmt_ci <- function(low, high) {
  paste0("[", fmt_num(low, 2), ", ", fmt_num(high, 2), "]")
}

fmt_est_se <- function(est, se, p = NULL) {
  if (is.null(p)) {
    paste0(fmt_num(est, 2), " (", fmt_num(se, 2), ")")
  } else {
    paste0(fmt_num(est, 2), stars_from_p(p), " (", fmt_num(se, 2), ")")
  }
}

show_table <- function(x, table_name) {
  cat("\n\n", table_name, "\n", sep = "")
  print(x, n = Inf, width = Inf)
  if (isTRUE(USE_VIEWER)) {
    View(x, title = table_name)
  }
  invisible(x)
}

make_reg_formula <- function(outcome, controls = control_vars) {
  as.formula(
    paste(
      outcome,
      "~ book_f +",
      paste(controls, collapse = " + ")
    )
  )
}

make_interaction_formula <- function(outcome) {
  as.formula(
    paste(
      outcome,
      "~ book_f * parent_edu_f +",
      paste(control_vars_no_parent_edu, collapse = " + ")
    )
  )
}


# ==============================================================================
# 4. Rubin 결합규칙 함수
# ==============================================================================

pool_svyglm_models <- function(model_list, df_com, outcome_name = "문해력", n_used) {

  coef_names <- Reduce(intersect, lapply(model_list, function(x) names(coef(x))))

  if (length(coef_names) == 0) {
    stop("결합할 공통 회귀계수가 없습니다. 범주 수준 또는 모형식을 확인하세요.")
  }

  coef_mat <- do.call(
    rbind,
    lapply(model_list, function(x) coef(x)[coef_names])
  )
  colnames(coef_mat) <- coef_names

  vcov_list <- lapply(model_list, function(x) {
    as.matrix(vcov(x))[coef_names, coef_names, drop = FALSE]
  })

  k <- length(model_list)
  q_bar <- colMeans(coef_mat)
  u_bar <- Reduce("+", vcov_list) / k

  if (k > 1) {
    b_mat <- stats::cov(coef_mat)
    if (is.null(dim(b_mat))) {
      b_mat <- matrix(b_mat, nrow = 1, dimnames = list(coef_names, coef_names))
    }
  } else {
    b_mat <- matrix(0, nrow = length(q_bar), ncol = length(q_bar), dimnames = list(coef_names, coef_names))
  }

  total_var_mat <- u_bar + (1 + 1 / k) * b_mat
  std_error <- sqrt(diag(total_var_mat))

  b_diag <- diag(b_mat)
  t_diag <- diag(total_var_mat)
  lambda <- ((1 + 1 / k) * b_diag) / t_diag
  lambda <- pmin(pmax(lambda, 0), 0.999999)

  df_old <- (k - 1) / (lambda^2)
  df_obs <- ((df_com + 1) / (df_com + 3)) * df_com * (1 - lambda)
  df_pooled <- 1 / ((1 / df_old) + (1 / df_obs))
  df_pooled[is.na(df_pooled) | is.infinite(df_pooled)] <- df_com
  df_pooled[b_diag < 1e-12] <- df_com

  statistic <- q_bar / std_error
  p_value <- 2 * pt(-abs(statistic), df = df_pooled)
  conf_low <- q_bar - qt(0.975, df = df_pooled) * std_error
  conf_high <- q_bar + qt(0.975, df = df_pooled) * std_error

  tibble(
    outcome = outcome_name,
    n = n_used,
    term = names(q_bar),
    estimate = as.numeric(q_bar),
    std_error = as.numeric(std_error),
    conf_low = as.numeric(conf_low),
    conf_high = as.numeric(conf_high),
    statistic_type = "t",
    statistic = as.numeric(statistic),
    df = as.numeric(df_pooled),
    p_value = as.numeric(p_value),
    stars = stars_from_p(p_value)
  )
}

pool_scalar_df <- function(df,
                           key_cols,
                           estimate_col,
                           se_col,
                           df_com,
                           outcome_name = "문해력",
                           n_used) {

  df_valid <- df %>%
    filter(!is.na(.data[[estimate_col]]), !is.na(.data[[se_col]]))

  if (length(key_cols) == 0) {
    grouped <- df_valid %>% mutate(.pool_key = "all") %>% group_by(.pool_key)
    select_keys <- character(0)
  } else {
    grouped <- df_valid %>% group_by(across(all_of(key_cols)))
    select_keys <- key_cols
  }

  pooled <- grouped %>%
    summarise(
      k = n(),
      estimate = mean(.data[[estimate_col]]),
      u_bar = mean((.data[[se_col]])^2),
      b_var = if (n() > 1) stats::var(.data[[estimate_col]]) else 0,
      .groups = "drop"
    )

  if (length(key_cols) == 0) {
    pooled <- pooled %>% select(-.pool_key)
  }

  pooled %>%
    mutate(
      total_var = u_bar + (1 + 1 / k) * b_var,
      std_error = sqrt(total_var),
      lambda = ((1 + 1 / k) * b_var) / total_var,
      lambda = pmin(pmax(lambda, 0), 0.999999),
      df_old = (k - 1) / (lambda^2),
      df_obs = ((df_com + 1) / (df_com + 3)) * df_com * (1 - lambda),
      df = 1 / ((1 / df_old) + (1 / df_obs)),
      df = if_else(is.na(df) | is.infinite(df) | b_var < 1e-12, df_com, df),
      conf_low = estimate - qt(0.975, df = df) * std_error,
      conf_high = estimate + qt(0.975, df = df) * std_error,
      statistic_type = "t",
      statistic = estimate / std_error,
      p_value = 2 * pt(-abs(statistic), df = df),
      stars = stars_from_p(p_value),
      outcome = outcome_name,
      n = n_used
    ) %>%
    select(
      outcome,
      n,
      all_of(select_keys),
      estimate,
      std_error,
      conf_low,
      conf_high,
      statistic_type,
      statistic,
      df,
      p_value,
      stars
    )
}


# ==============================================================================
# 5. 데이터 불러오기 및 전처리
# ==============================================================================

if (!file.exists(DATA_PATH)) {
  stop("DATA_PATH에 지정한 파일이 없습니다: ", DATA_PATH)
}

dat_raw <- readr::read_delim(
  DATA_PATH,
  delim = NULL,
  show_col_types = FALSE,
  locale = readr::locale(encoding = "UTF-8")
)

# 원자료에서 필요한 변수가 있는지 확인
required_raw_vars <- c(
  pv_lit,
  main_weight,
  rep_weights,
  book_var,
  age_var,
  "GENDER_R",
  "PAREDC2",
  "J2_Q04d",
  "J2_Q05d",
  "J2_Q07_C",
  "J2_Q0801",
  "J2_Q0802"
)

missing_raw_vars <- setdiff(required_raw_vars, names(dat_raw))
if (length(missing_raw_vars) > 0) {
  stop("원자료에서 필요한 변수를 찾지 못했습니다: ", paste(missing_raw_vars, collapse = ", "))
}

# 문자형 특수결측은 전체 문자 변수에서 처리
# 숫자형 특수결측은 원문 설문변수 중심으로 처리
special_numeric_vars <- c(
  book_var,
  "PAREDC2",
  "J2_Q04d",
  "J2_Q05d",
  "J2_Q07_C",
  "J2_Q0801",
  "J2_Q0802"
)

dat <- dat_raw %>%
  mutate(across(where(is.character), recode_char_missing)) %>%
  mutate(across(all_of(special_numeric_vars), recode_numeric_special_missing)) %>%
  mutate(
    # 도서 보유 수
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

    # 연령 범주
    age_f = case_when(
      as_code(.data[[age_var]]) == 1 ~ "20-24세",
      as_code(.data[[age_var]]) == 2 ~ "25-29세",
      as_code(.data[[age_var]]) == 3 ~ "30-34세",
      as_code(.data[[age_var]]) == 4 ~ "35-39세",
      as_code(.data[[age_var]]) == 5 ~ "40-44세",
      as_code(.data[[age_var]]) == 6 ~ "45-49세",
      as_code(.data[[age_var]]) == 7 ~ "50-54세",
      as_code(.data[[age_var]]) == 8 ~ "55-59세",
      as_code(.data[[age_var]]) == 9 ~ "60-64세",
      as_code(.data[[age_var]]) == 10 ~ "65세 이상",
      TRUE ~ NA_character_
    ),
    age_f = factor(
      age_f,
      levels = c(
        "20-24세", "25-29세", "30-34세", "35-39세",
        "40-44세", "45-49세", "50-54세", "55-59세",
        "60-64세", "65세 이상"
      )
    ),

    # 부모 교육수준
    parent_edu_f = case_when(
      as_code(PAREDC2) == 1 ~ "중졸 이하",
      as_code(PAREDC2) == 2 ~ "고졸/전문대졸",
      as_code(PAREDC2) == 3 ~ "대졸 이상",
      TRUE ~ NA_character_
    ),
    parent_edu_f = factor(parent_edu_f, levels = c("중졸 이하", "고졸/전문대졸", "대졸 이상")),

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
    residence_f = factor(residence_f, levels = c("대도시", "중소도시", "소도시/읍면", "농어촌/시골")),

    # 14세 당시 가족구조
    family14_f = case_when(
      as_code(J2_Q0801) == 1 & as_code(J2_Q0802) == 1 ~ "양부모 동거",
      as_code(J2_Q0801) == 1 & as_code(J2_Q0802) != 1 ~ "모만 동거",
      as_code(J2_Q0801) != 1 & as_code(J2_Q0802) == 1 ~ "부만 동거",
      as_code(J2_Q0801) != 1 & as_code(J2_Q0802) != 1 ~ "생부모 비동거",
      TRUE ~ NA_character_
    ),
    family14_f = factor(family14_f, levels = c("양부모 동거", "모만 동거", "부만 동거", "생부모 비동거"))
  )

# 기준범주 설정
dat <- dat %>%
  mutate(
    book_f = relevel(book_f, ref = "10권 이하"),
    parent_edu_f = relevel(parent_edu_f, ref = "중졸 이하")
  )

# 완전사례분석에 필요한 변수
analysis_vars <- c(
  pv_lit,
  "book_f",
  control_vars,
  main_weight,
  rep_weights
)

# 완전사례 자료 생성
cc_data <- dat %>%
  select(all_of(analysis_vars)) %>%
  tidyr::drop_na() %>%
  mutate(across(where(is.factor), droplevels)) %>%
  mutate(
    book_f = relevel(factor(book_f, levels = book_levels), ref = "10권 이하"),
    parent_edu_f = relevel(factor(parent_edu_f, levels = c("중졸 이하", "고졸/전문대졸", "대졸 이상")), ref = "중졸 이하")
  )

cc_n <- nrow(cc_data)
cc_des <- make_design(cc_data)
cc_df <- survey::degf(cc_des)

cat("\n완전사례분석 표본 크기: ", cc_n, "명\n", sep = "")
cat("복합표본 설계 자유도: ", cc_df, "\n", sep = "")


# ==============================================================================
# 6. 표 1. 도서 보유 수별 외생적 배경변수 기술통계 및 교차분석
# ==============================================================================

background_specs <- list(
  list(
    var = "gender_f",
    label = "성별",
    levels = c("남성", "여성")
  ),
  list(
    var = "age_f",
    label = "연령",
    levels = c(
      "20-24세",
      "25-29세",
      "30-34세",
      "35-39세",
      "40-44세",
      "45-49세",
      "50-54세",
      "55-59세",
      "60-64세",
      "65세 이상"
    )
  ),
  list(
    var = "parent_edu_f",
    label = "부모 교육수준",
    levels = c("중졸 이하", "고졸/전문대졸", "대졸 이상")
  ),
  list(
    var = "mother_work_f",
    label = "모 경제활동",
    levels = c("유급직", "무직/가사 등")
  ),
  list(
    var = "father_work_f",
    label = "부 경제활동",
    levels = c("유급직", "무직/가사 등")
  ),
  list(
    var = "residence_f",
    label = "거주지",
    levels = c("대도시", "중소도시", "소도시/읍면", "농어촌/시골")
  ),
  list(
    var = "family14_f",
    label = "가족구조",
    levels = c("양부모 동거", "모만 동거", "부만 동거", "생부모 비동거")
  )
)

safe_svy_chisq_p <- function(var) {
  fml <- as.formula(paste("~ book_f +", var))

  try_methods <- c("F", "adjWald", "Chisq")

  for (method in try_methods) {
    out <- tryCatch(
      survey::svychisq(fml, cc_des, statistic = method),
      error = function(e) NULL,
      warning = function(w) suppressWarnings(survey::svychisq(fml, cc_des, statistic = method))
    )

    if (!is.null(out)) {
      p <- suppressWarnings(as.numeric(out$p.value))
      if (length(p) == 1 && is.finite(p)) {
        return(p)
      }
    }
  }

  NA_real_
}

weighted_col_percent <- function(var, category) {
  tmp <- cc_data %>%
    filter(!is.na(book_f), !is.na(.data[[var]])) %>%
    mutate(category_tmp = as.character(.data[[var]])) %>%
    group_by(book_f, category_tmp) %>%
    summarise(weighted_n = sum(.data[[main_weight]], na.rm = TRUE), .groups = "drop") %>%
    group_by(book_f) %>%
    mutate(percent = 100 * weighted_n / sum(weighted_n, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(category_tmp == category) %>%
    select(book_f, percent)

  out <- setNames(rep(NA_real_, length(book_levels)), book_levels)
  out[as.character(tmp$book_f)] <- tmp$percent
  out
}

make_table1_section <- function(spec) {
  p_val <- safe_svy_chisq_p(spec$var)

  header <- tibble(
    변수명 = spec$label,
    `≤10` = "",
    `11-25` = "",
    `26-100` = "",
    `101-200` = "",
    `201-500` = "",
    `>500` = "",
    `pᵃ` = fmt_p_with_stars(p_val)
  )

  rows <- map_dfr(spec$levels, function(cat_label) {
    pct <- weighted_col_percent(spec$var, cat_label)
    tibble(
      변수명 = cat_label,
      `≤10` = fmt_num(pct["10권 이하"], 1),
      `11-25` = fmt_num(pct["11-25권"], 1),
      `26-100` = fmt_num(pct["26-100권"], 1),
      `101-200` = fmt_num(pct["101-200권"], 1),
      `201-500` = fmt_num(pct["201-500권"], 1),
      `>500` = fmt_num(pct["500권 초과"], 1),
      `pᵃ` = ""
    )
  })

  bind_rows(header, rows)
}

표1_기술통계_교차분석 <- map_dfr(background_specs, make_table1_section)

show_table(
  표1_기술통계_교차분석,
  "<표 1> 14세 시점 도서 보유 수별 외생적 배경변수의 기술통계 및 교차분석 결과"
)

cat("\n표 1 주. 수치는 복합표본 최종가중치를 적용한 열 백분율(%)임; a Rao-Scott 수정 카이제곱 검정의 p-value임.\n")


# ==============================================================================
# 7. 표 2. 도서 보유 수와 성인기 문해력의 회귀분석
# ==============================================================================

fit_pv_regression_models <- function(formula_fun) {
  map(pv_lit, function(y) {
    survey::svyglm(
      formula_fun(y),
      design = cc_des
    )
  })
}

main_model_list <- fit_pv_regression_models(function(y) make_reg_formula(y, control_vars))
main_pooled <- pool_svyglm_models(main_model_list, df_com = cc_df, n_used = cc_n)

main_effect_results <- main_pooled %>%
  filter(str_detect(term, "^book_f")) %>%
  mutate(도서보유수 = str_remove(term, "^book_f"))

표2_회귀분석 <- main_effect_results %>%
  transmute(
    `도서 보유 수` = 도서보유수,
    추정치 = fmt_num(estimate, 2),
    표준오차 = fmt_num(std_error, 2),
    t = fmt_num(statistic, 2),
    p = fmt_p_with_stars(p_value),
    `95% 신뢰구간` = fmt_ci(conf_low, conf_high)
  )

show_table(
  표2_회귀분석,
  "<표 2> 14세 시점 도서 보유 수와 성인기 문해력의 관계"
)

cat("\n표 2 주. 표본 크기는 ", cc_n, "명이다. 기준집단은 ‘10권 이하’이다. ", sep = "")
cat("모형은 성별, 연령, 부모 교육수준, 14세 당시 부모의 경제활동(유급직), 거주지역 규모, 가족구조를 통제하였다. ")
cat("모수 및 표준오차는 10개의 Plausible Values 각각에 대하여 복합표본 설계(최종가중치 및 80개의 반복가중치를 적용한 Fay의 BRR, 계수: 0.3)를 반영해 추정한 뒤, Rubin의 결합규칙에 따라 최종 통합되었다. ")
cat("† p < .10, * p < .05, ** p < .01, *** p < .001.\n")


# ==============================================================================
# 8. 표 3. 민감도 분석
# ==============================================================================

extract_sense_rv <- function(sense_out) {
  ss <- sense_out$sensitivity_stats

  rv <- NA_real_
  rv_alpha <- NA_real_

  if (!is.null(ss)) {
    ss_names <- names(ss)

    if ("rv_q" %in% ss_names) rv <- as.numeric(ss[["rv_q"]])
    if ("rv_qa" %in% ss_names) rv_alpha <- as.numeric(ss[["rv_qa"]])

    # 패키지 버전에 따른 대체 이름 보정
    if (is.na(rv) && "RV_q" %in% ss_names) rv <- as.numeric(ss[["RV_q"]])
    if (is.na(rv_alpha) && "RV_qa" %in% ss_names) rv_alpha <- as.numeric(ss[["RV_qa"]])
  }

  tibble(RV = rv, 극단적_RV = rv_alpha)
}

run_sensitivity_pv <- function(i) {
  fml <- make_reg_formula(pv_lit[i], control_vars)

  fit_svy <- survey::svyglm(fml, design = cc_des)
  fit_sum <- coef(summary(fit_svy))

  if (!(SENSE_TREAT_TERM %in% rownames(fit_sum))) {
    stop("민감도 분석 대상 계수가 모형에 없습니다: ", SENSE_TREAT_TERM)
  }

  est_brr <- fit_sum[SENSE_TREAT_TERM, "Estimate"]
  se_brr <- fit_sum[SENSE_TREAT_TERM, "Std. Error"]

  # sensemakr는 lm 객체를 요구하므로, 민감도 계산에는 최종가중치를 적용한 lm을 사용한다.
  # 조정 추정치와 조정 표준오차는 위에서 추정한 복합표본 추정치와 표준오차를 기준으로 다시 계산한다.
  fit_lm <- stats::lm(
    fml,
    data = cc_data,
    weights = cc_data[[main_weight]]
  )

  sense_out <- sensemakr::sensemakr(
    model = fit_lm,
    treatment = SENSE_TREAT_TERM,
    benchmark_covariates = SENSE_BENCHMARK_COV,
    kd = SENSE_KD,
    q = 1,
    alpha = 0.05
  )

  bounds <- as_tibble(sense_out$bounds) %>%
    mutate(
      scenario = case_when(
        str_detect(bound_label, "1") ~ "1배",
        str_detect(bound_label, "2") ~ "2배",
        str_detect(bound_label, "3") ~ "3배",
        TRUE ~ as.character(bound_label)
      ),
      Adjusted_Estimate = sensemakr::adjusted_estimate(
        estimate = est_brr,
        se = se_brr,
        dof = cc_df,
        r2dz.x = r2dz.x,
        r2yz.dx = r2yz.dx
      ),
      Adjusted_SE = sensemakr::adjusted_se(
        se = se_brr,
        dof = cc_df,
        r2dz.x = r2dz.x,
        r2yz.dx = r2yz.dx
      ),
      pv = i
    )

  rv <- extract_sense_rv(sense_out) %>% mutate(pv = i)

  list(bounds = bounds, rv = rv)
}

sensitivity_pv_results <- map(seq_along(pv_lit), run_sensitivity_pv)

sensitivity_bounds_all <- map_dfr(sensitivity_pv_results, "bounds")
sensitivity_rv_all <- map_dfr(sensitivity_pv_results, "rv")

sensitivity_bounds_pooled <- pool_scalar_df(
  df = sensitivity_bounds_all,
  key_cols = "scenario",
  estimate_col = "Adjusted_Estimate",
  se_col = "Adjusted_SE",
  df_com = cc_df,
  n_used = cc_n
)

sensitivity_rv_pooled <- sensitivity_rv_all %>%
  summarise(
    RV = mean(RV, na.rm = TRUE),
    극단적_RV = mean(극단적_RV, na.rm = TRUE)
  )

표3_민감도분석 <- bind_rows(
  tibble(
    구분 = "강건성 지표",
    RV = fmt_num(sensitivity_rv_pooled$RV, 3),
    `극단적 RV` = fmt_num(sensitivity_rv_pooled$극단적_RV, 3),
    `조정 추정치(조정 표준오차)` = "—"
  ),
  sensitivity_bounds_pooled %>%
    transmute(
      구분 = scenario,
      RV = "—",
      `극단적 RV` = "—",
      `조정 추정치(조정 표준오차)` = paste0(fmt_num(estimate, 2), "(", fmt_num(std_error, 2), ")")
    )
)

show_table(
  표3_민감도분석,
  "<표 3> 미관측 교란요인에 대한 도서 보유 수 효과의 민감도 분석"
)

cat("\n표 3 주. 표본 크기는 ", cc_n, "명이다. 기준집단은 ‘10권 이하’이며, 분석 대상 효과는 ‘201-500권’ 집단의 효과이다. ", sep = "")
cat("RV는 미관측 교란요인이 현재 추정 효과를 통계적으로 유의하지 않게 만들기 위해 필요한 최소 설명력을 의미한다. ")
cat("극단적 RV는 결과변수와의 관련성을 제한하지 않는 극단적 시나리오에서, 미관측 교란요인이 원인변수와 가져야 하는 최소 설명력을 의미한다. ")
cat("1배, 2배, 3배 시나리오는 미관측 교란요인이 부모 교육수준(대졸 이상) 대비 각각 1배, 2배, 3배 강력한 설명력을 가진다고 가정한 경우이다. ")
cat("조정치는 10개의 PVs 각각에 복합표본 설계(Fay의 BRR, 계수: 0.3)를 반영한 뒤 Rubin의 결합규칙에 따라 통합하였다.\n")


# ==============================================================================
# 9. 표 4. 도서 보유 수와 부모 교육수준의 상호작용항
# ==============================================================================

interaction_model_list <- fit_pv_regression_models(function(y) make_interaction_formula(y))
interaction_pooled <- pool_svyglm_models(interaction_model_list, df_com = cc_df, n_used = cc_n)

interaction_results <- interaction_pooled %>%
  filter(str_detect(term, ":")) %>%
  mutate(
    항 = term %>%
      str_replace("book_f", "도서 ") %>%
      str_replace("parent_edu_f", "부모교육 ") %>%
      str_replace(":", " × ")
  )

표4_상호작용항 <- interaction_results %>%
  transmute(
    상호작용항 = 항,
    추정치 = fmt_num(estimate, 2),
    표준오차 = fmt_num(std_error, 2),
    t = fmt_num(statistic, 2),
    p = fmt_p_with_stars(p_value),
    `95% 신뢰구간` = fmt_ci(conf_low, conf_high)
  )

show_table(
  표4_상호작용항,
  "<표 4> 도서 보유 수와 부모 교육수준의 상호작용항 분석"
)

cat("\n표 4 주. 표본 크기는 ", cc_n, "명이다. 상호작용항의 기준 범주는 ‘부모 교육수준 중졸 이하 × 도서 10권 이하’이다. ", sep = "")
cat("모형은 도서 보유 수, 부모 교육수준, 도서 보유 수와 부모 교육수준의 상호작용항을 포함하며, 성별, 연령, 14세 당시 부모의 경제활동 여부, 거주지역 규모, 가족구조를 통제하였다. ")
cat("추정치 및 표준오차는 10개의 PVs 각각에 대해 복합표본 설계, 즉 최종가중치와 80개의 반복가중치를 적용한 Fay의 BRR, 계수 0.3을 반영하여 산출한 뒤 Rubin의 결합규칙에 따라 통합하였다. † p < .10, * p < .05, ** p < .01, *** p < .001.\n")


# ==============================================================================
# 10. 표 5. 도서풍부 취약계층과 도서빈약 상위배경 집단의 사후 대비
# ==============================================================================

# 표 5는 도서보유수 × 부모교육수준 상호작용 모형에서 산출한
# 추정주변평균(estimated marginal means) 간 선택 대비 분석이다.
# 값은 기준집단 - 비교집단이다.
# 각 PV별로 복합표본 설계(Fay's BRR)를 반영한 뒤 Rubin 결합규칙으로 통합한다.

표5_대비설정 <- tibble::tribble(
  ~contrast_id, ~기준_부모교육수준, ~기준_도서보유수, ~비교_부모교육수준, ~비교_도서보유수,
  
  "C1",
  "중졸 이하", "201-500권",
  "고졸/전문대졸", "10권 이하",
  
  "C2",
  "중졸 이하", "201-500권",
  "고졸/전문대졸", "11-25권",
  
  "C3",
  "중졸 이하", "201-500권",
  "대졸 이상", "10권 이하",
  
  "C4",
  "중졸 이하", "201-500권",
  "대졸 이상", "11-25권"
) %>%
  mutate(
    기준집단 = paste0("부모 ", 기준_부모교육수준, " & 도서 ", 기준_도서보유수),
    비교집단 = paste0("부모 ", 비교_부모교육수준, " & 도서 ", 비교_도서보유수),
    contrast_id = factor(contrast_id, levels = c("C1", "C2", "C3", "C4"))
  )

get_table5_one_pv <- function(pv_name, pv_index) {
  
  # 기존 코드에서 정의한 상호작용 모형 사용:
  # PVLIT ~ book_f * parent_edu_f + 성별 + 연령 + 부모경제활동 + 거주지 + 가족구조
  fml <- make_interaction_formula(pv_name)
  
  fit <- survey::svyglm(
    fml,
    design = cc_des
  )
  
  fit$data <- cc_des$variables
  
  emm <- emmeans::emmeans(
    fit,
    ~ book_f * parent_edu_f,
    data = cc_des$variables
  )
  
  emm_df <- as.data.frame(emm) %>%
    as_tibble() %>%
    mutate(
      도서보유수 = as.character(book_f),
      부모교육수준 = as.character(parent_edu_f)
    )
  
  purrr::map_dfr(seq_len(nrow(표5_대비설정)), function(r) {
    
    row <- 표5_대비설정[r, ]
    
    idx_base <- which(
      emm_df$부모교육수준 == row$기준_부모교육수준 &
        emm_df$도서보유수 == row$기준_도서보유수
    )
    
    idx_comp <- which(
      emm_df$부모교육수준 == row$비교_부모교육수준 &
        emm_df$도서보유수 == row$비교_도서보유수
    )
    
    if (length(idx_base) != 1 | length(idx_comp) != 1) {
      cat("\n현재 emmeans 셀 목록:\n")
      print(
        emm_df %>%
          select(도서보유수, 부모교육수준),
        n = Inf
      )
      
      stop(
        "선택 대비에 필요한 집단을 찾지 못했습니다: ",
        row$contrast_id
      )
    }
    
    # 기준집단 - 비교집단
    contrast_vector <- rep(0, nrow(emm_df))
    contrast_vector[idx_base] <- 1
    contrast_vector[idx_comp] <- -1
    
    contr <- emmeans::contrast(
      emm,
      method = list(선택대비 = contrast_vector),
      adjust = "none"
    )
    
    as.data.frame(contr) %>%
      as_tibble() %>%
      transmute(
        contrast_id = row$contrast_id,
        기준집단 = row$기준집단,
        비교집단 = row$비교집단,
        pv = pv_index,
        estimate = estimate,
        SE = SE
      )
  })
}

표5_raw <- purrr::map2_dfr(
  pv_lit,
  seq_along(pv_lit),
  ~ get_table5_one_pv(.x, .y)
)

pool_table5_contrasts <- function(df, df_com, n_used) {
  
  df %>%
    mutate(
      estimate = as.numeric(estimate),
      SE = as.numeric(SE)
    ) %>%
    filter(
      !is.na(estimate),
      !is.na(SE),
      is.finite(estimate),
      is.finite(SE)
    ) %>%
    group_by(contrast_id, 기준집단, 비교집단) %>%
    summarise(
      k = n(),
      estimate = mean(estimate, na.rm = TRUE),
      u_bar = mean(SE^2, na.rm = TRUE),
      b_var = if_else(
        k > 1,
        stats::var(estimate, na.rm = TRUE),
        0
      ),
      .groups = "drop"
    ) %>%
    mutate(
      b_var = if_else(is.na(b_var), 0, b_var),
      total_var = u_bar + (1 + 1 / k) * b_var,
      std_error = sqrt(total_var),
      
      lambda = ((1 + 1 / k) * b_var) / total_var,
      lambda = if_else(is.na(lambda) | is.infinite(lambda), 0, lambda),
      lambda = pmin(pmax(lambda, 0), 0.999999),
      
      df_old = if_else(lambda > 0, (k - 1) / lambda^2, Inf),
      df_obs = ((df_com + 1) / (df_com + 3)) * df_com * (1 - lambda),
      df = 1 / ((1 / df_old) + (1 / df_obs)),
      df = if_else(is.na(df) | is.infinite(df) | b_var < 1e-12, df_com, df),
      
      statistic = estimate / std_error,
      p_value = 2 * pt(-abs(statistic), df = df),
      conf_low = estimate - qt(0.975, df = df) * std_error,
      conf_high = estimate + qt(0.975, df = df) * std_error,
      
      p_adjust_method = POSTHOC_ADJUST_METHOD,
      p_adjusted = p.adjust(p_value, method = POSTHOC_ADJUST_METHOD),
      stars = stars_from_p(p_value),
      stars_adjusted = stars_from_p(p_adjusted),
      n = n_used
    )
}

표5_사후대비_수치 <- pool_table5_contrasts(
  df = 표5_raw,
  df_com = cc_df,
  n_used = cc_n
) %>%
  arrange(contrast_id)

표5_사후대비 <- 표5_사후대비_수치 %>%
  transmute(
    기준집단,
    비교집단,
    `추정주변평균 차이` = fmt_num(estimate, 2),
    표준오차 = fmt_num(std_error, 2),
    t = fmt_num(statistic, 2),
    p = p_format(p_value),
    `보정 p` = p_format(p_adjusted),
    유의성 = stars_adjusted,
    `95% 신뢰구간` = fmt_ci(conf_low, conf_high)
  )

show_table(
  표5_사후대비,
  "<표 5> 문해력에서 도서풍부 취약계층과 도서빈약 상위배경 집단의 사후 대비 결과"
)

cat("\n표 5 주. 표의 값은 기준집단에서 비교집단을 뺀 추정주변평균 차이이다. ")
cat("양수 값은 기준집단의 문해력 점수가 비교집단보다 높음을 의미한다. ")
cat("모형은 도서 보유 수와 부모 교육수준의 상호작용항을 포함하며, 성별, 연령, 14세 당시 부모의 경제활동 여부, 거주지역 규모, 가족구조를 통제하였다. ")
cat("각 대비는 10개의 Plausible Values 각각에 대해 복합표본 설계, 즉 최종가중치와 80개의 반복가중치를 적용한 Fay의 BRR, 계수 0.3을 반영하여 산출한 뒤 Rubin의 결합규칙으로 통합하였다. ")
cat("보정 p값은 ", POSTHOC_ADJUST_METHOD, " 방식으로 산출하였다. † p < .10, * p < .05, ** p < .01, *** p < .001.\n", sep = "")


################################################################################
# 끝
################################################################################
