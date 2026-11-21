# PIAAC 2주기 한국 자료: 청소년기 도서 보유 수와 성인기 문해력
# 실행: 작업 폴더에 prgkorp2.csv를 두고 이 파일 전체를 실행합니다.
# 출력: output_piaac_sensitivity 폴더에 표 1~3, 상세 결과 및 실행 환경을 저장합니다.
# 계산 가정과 함수별 처리 절차는 해당 코드의 주석을 참조하세요.


# 0. 실행 옵션
# USE_VIEWER: 대화형 실행 시 표를 별도 창에 표시합니다.
# SAVE_OUTPUT: CSV, RDS 및 session_info.txt의 저장 여부를 지정합니다.
DATA_PATH   <- "prgkorp2.csv"
USE_VIEWER  <- interactive()
SAVE_OUTPUT <- TRUE
OUTPUT_DIR  <- "output_piaac_sensitivity"


# 1. 패키지: 설치되지 않은 패키지는 설치한 뒤 불러옵니다.
packages <- c("tidyverse", "survey", "sensemakr")

to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

invisible(lapply(packages, library, character.only = TRUE))


# 2. 분석 설정: 문해력 PV 10개, 최종가중치, Fay 반복가중치 80개를 사용합니다.
pv_lit      <- paste0("PVLIT", 1:10)
main_weight <- "SPFWT0"
rep_weights <- paste0("SPFWT", 1:80)
fay_rho     <- 0.3

book_var <- "J2_Q06"
age_var  <- "AGEG5LFS"

book_levels <- c(
  "10권 이하", "11-25권", "26-100권",
  "101-200권", "201-500권", "500권 초과"
)
book_ref <- "10권 이하"

book_col_labels <- c(
  "10권 이하"  = "≤10",
  "11-25권"    = "11-25",
  "26-100권"   = "26-100",
  "101-200권"  = "101-200",
  "201-500권"  = "201-500",
  "500권 초과" = ">500"
)

control_vars <- c(
  "gender_f",
  "age_f",
  "parent_edu_f",
  "mother_work_f",
  "father_work_f",
  "residence_f",
  "family14_f"
)

# 각 비기준 도서 범주의 회귀계수를 개별적으로 평가합니다.
# SENSE_KD는 비교기준 보정식의 원인·결과 양쪽 배수(kd, ky)에 공통 적용합니다.
# SENSE_ALPHA는 민감도 분석의 양측 유의수준입니다.
SENSE_TREAT_TERMS   <- paste0("book_f", setdiff(book_levels, book_ref))
SENSE_BENCHMARK_COV <- "parent_edu_f대졸 이상"
SENSE_KD            <- c(1, 2, 3)
SENSE_ALPHA         <- 0.05

scenario_levels <- paste0(SENSE_KD, "배")


# 3. 보조 함수: 범주 코드를 정수로 변환합니다.
as_code <- function(x) suppressWarnings(as.integer(as.character(x)))

# 문자형 특수 결측 표기를 NA로 통일합니다.
recode_char_missing <- function(x) {
  if (is.character(x)) x[tolower(trimws(x)) %in% c(".", ".n", ".r", ".d", ".v", ".a", ".b", "", "na")] <- NA_character_
  x
}

# 지정한 범주형 변수에서만 음수 및 특수 결측 코드를 NA로 바꿉니다.
recode_numeric_special_missing <- function(x) {
  x_num <- as_code(x)
  x_num[x_num < 0 | x_num %in% c(7, 8, 9, 96, 97, 98, 99, 996, 997, 998, 999)] <- NA_integer_
  x_num
}

# 최종가중치와 제공된 반복가중치로 Fay BRR 설계 객체를 만듭니다.
# combined.weights=TRUE: 반복가중치는 최종가중치를 포함한 완전한 가중치입니다.
# mse=TRUE: 반복분산을 반복추정치와 최종가중치 추정치의 차이로 계산합니다.
make_design <- function(data) {
  survey::svrepdesign(
    weights          = as.formula(paste0("~", main_weight)),
    repweights       = data[, rep_weights],
    data             = data,
    type             = "Fay",
    rho              = fay_rho,
    combined.weights = TRUE,
    mse              = TRUE
  )
}

# PV 이름과 통제변수 목록으로 동일한 회귀모형식을 구성합니다.
make_reg_formula <- function(outcome, controls = control_vars) {
  as.formula(paste(outcome, "~ book_f +", paste(controls, collapse = " + ")))
}

# 표 출력용 함수: 소수 자릿수, p값, 별표 및 신뢰구간 표기를 정리합니다.
fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

fmt_dec <- function(x, digits = 3) {
  ifelse(is.na(x), "", sub("^0", "", sprintf(paste0("%.", digits, "f"), x)))
}

p_format <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "<.001",
    TRUE      ~ fmt_dec(p, 3)
  )
}

stars_from_p <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

fmt_p_with_stars <- function(p) paste0(p_format(p), stars_from_p(p))

fmt_ci <- function(low, high) paste0("[", fmt_num(low, 2), ", ", fmt_num(high, 2), "]")

# 표를 콘솔에 출력하고, 실행 옵션에 따라 뷰어 표시 및 CSV 저장을 수행합니다.
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


# 4. PV별 회귀모형의 추정치와 공분산을 Rubin의 결합규칙으로 통합합니다.
# 총분산 = PV 내 분산의 평균 + (1 + 1/k) × PV 간 추정치 분산.
# df_com은 survey::degf로 구한 설계 자유도이며, Barnard–Rubin 형태로 보정합니다.
# PV 간 분산이 사실상 0이면 df_com을 사용하고, 양측 t 검정과 95% CI를 계산합니다.
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
  q_bar <- colMeans(coef_mat)
  u_bar <- Reduce("+", vcov_list) / k

  if (k > 1) {
    b_mat <- stats::cov(coef_mat)
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
  lambda <- pmin(pmax(((1 + 1 / k) * b_diag) / diag(total_var_mat), 0), 0.999999)

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

# PV별 단일 추정치와 표준오차를 key_cols별로 통합합니다.
# 회귀계수 통합과 동일한 Rubin 분산 및 자유도 보정을 조정 추정치에 적용합니다.
pool_scalar_df <- function(df, key_cols, estimate_col, se_col,
                           df_com, outcome_name = "문해력", n_used) {

  stopifnot(length(key_cols) >= 1)

  df %>%
    rename(.est = !!as.name(estimate_col), .se = !!as.name(se_col)) %>%
    filter(!is.na(.est), !is.na(.se)) %>%
    group_by(across(all_of(key_cols))) %>%
    summarise(
      k        = n(),
      estimate = mean(.est),
      u_bar    = mean(.se^2),
      b_var    = if (n() > 1) stats::var(.est) else 0,
      .groups  = "drop"
    ) %>%
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
      outcome, n, all_of(key_cols),
      estimate, std_error, conf_low, conf_high,
      statistic_type, statistic, df, p_value, stars
    )
}


# 5. 자료 읽기: 첫 줄에서 가장 많이 나타나는 구분자를 선택합니다.
# 지원 구분자는 세미콜론, 쉼표, 탭, 세로줄이며 UTF-8로 읽습니다.
read_piaac <- function(path) {
  if (!file.exists(path)) {
    stop("DATA_PATH에 지정한 파일이 없습니다: ", path)
  }

  first_line <- readLines(path, n = 1, warn = FALSE)
  cands      <- c(";", ",", "\t", "|")
  counts     <- vapply(cands, function(d) {
    length(gregexpr(d, first_line, fixed = TRUE)[[1]][
      gregexpr(d, first_line, fixed = TRUE)[[1]] > 0])
  }, integer(1))

  if (max(counts) == 0) {
    stop("구분자를 판별하지 못했습니다. 파일 첫 줄을 확인하세요.")
  }
  delim <- cands[which.max(counts)]

  cat(sprintf("[자료] 구분자 '%s'로 읽습니다.\n",
              if (delim == "\t") "\\t" else delim))

  readr::read_delim(
    path,
    delim          = delim,
    show_col_types = FALSE,
    progress       = FALSE,
    locale         = readr::locale(encoding = "UTF-8")
  )
}

dat_raw <- read_piaac(DATA_PATH)

cat(sprintf("[자료] %s행 %s열\n",
            format(nrow(dat_raw), big.mark = ","),
            format(ncol(dat_raw), big.mark = ",")))

# 분석에 필요한 원변수가 모두 있는지 전처리 전에 확인합니다.
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

# 범주형 변수에만 특수 숫자 결측 처리를 적용합니다(PV와 가중치에는 미적용).
special_numeric_vars <- c(
  book_var, "PAREDC2", "J2_Q04d", "J2_Q05d",
  "J2_Q07_C", "J2_Q0801", "J2_Q0802"
)

# 원자료 코드를 분석용 범주로 변환합니다. factor 수준의 첫 범주가 기본 기준집단입니다.
dat <- dat_raw %>%
  mutate(across(where(is.character), recode_char_missing)) %>%
  mutate(across(all_of(special_numeric_vars), recode_numeric_special_missing)) %>%
  mutate(across(all_of(c(pv_lit, main_weight, rep_weights)),
                ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(
    book_f = factor(
      case_when(
        as_code(.data[[book_var]]) == 1 ~ "10권 이하",
        as_code(.data[[book_var]]) == 2 ~ "11-25권",
        as_code(.data[[book_var]]) == 3 ~ "26-100권",
        as_code(.data[[book_var]]) == 4 ~ "101-200권",
        as_code(.data[[book_var]]) == 5 ~ "201-500권",
        as_code(.data[[book_var]]) == 6 ~ "500권 초과",
        TRUE ~ NA_character_
      ),
      levels = book_levels
    ),

    gender_f = factor(
      case_when(
        as_code(GENDER_R) == 1 ~ "남성",
        as_code(GENDER_R) == 2 ~ "여성",
        TRUE ~ NA_character_
      ),
      levels = c("남성", "여성")
    ),

    age_f = factor(
      case_when(
        as_code(.data[[age_var]]) == 1  ~ "16-19세",
        as_code(.data[[age_var]]) == 2  ~ "20-24세",
        as_code(.data[[age_var]]) == 3  ~ "25-29세",
        as_code(.data[[age_var]]) == 4  ~ "30-34세",
        as_code(.data[[age_var]]) == 5  ~ "35-39세",
        as_code(.data[[age_var]]) == 6  ~ "40-44세",
        as_code(.data[[age_var]]) == 7  ~ "45-49세",
        as_code(.data[[age_var]]) == 8  ~ "50-54세",
        as_code(.data[[age_var]]) == 9  ~ "55-59세",
        as_code(.data[[age_var]]) == 10 ~ "60-65세",
        TRUE ~ NA_character_
      ),
      levels = c("16-19세", "20-24세", "25-29세", "30-34세", "35-39세",
                 "40-44세", "45-49세", "50-54세", "55-59세", "60-65세")
    ),

    parent_edu_f = factor(
      case_when(
        as_code(PAREDC2) == 1 ~ "중졸 이하",
        as_code(PAREDC2) == 2 ~ "고졸 수준",
        as_code(PAREDC2) == 3 ~ "대졸 이상",
        TRUE ~ NA_character_
      ),
      levels = c("중졸 이하", "고졸 수준", "대졸 이상")
    ),

    mother_work_f = factor(
      case_when(
        as_code(J2_Q04d) == 1 ~ "유급직",
        as_code(J2_Q04d) == 2 ~ "무직/가사 등",
        TRUE ~ NA_character_
      ),
      levels = c("유급직", "무직/가사 등")
    ),

    father_work_f = factor(
      case_when(
        as_code(J2_Q05d) == 1 ~ "유급직",
        as_code(J2_Q05d) == 2 ~ "무직/가사 등",
        TRUE ~ NA_character_
      ),
      levels = c("유급직", "무직/가사 등")
    ),

    residence_f = factor(
      case_when(
        as_code(J2_Q07_C) == 1 ~ "대도시",
        as_code(J2_Q07_C) == 2 ~ "중소도시",
        as_code(J2_Q07_C) == 3 ~ "소도시/읍면",
        as_code(J2_Q07_C) == 4 ~ "농어촌/시골",
        TRUE ~ NA_character_
      ),
      levels = c("대도시", "중소도시", "소도시/읍면", "농어촌/시골")
    ),

    family14_f = factor(
      case_when(
        as_code(J2_Q0801) == 1 & as_code(J2_Q0802) == 1 ~ "양부모 동거",
        as_code(J2_Q0801) == 1 & as_code(J2_Q0802) != 1 ~ "모만 동거",
        as_code(J2_Q0801) != 1 & as_code(J2_Q0802) == 1 ~ "부만 동거",
        as_code(J2_Q0801) != 1 & as_code(J2_Q0802) != 1 ~ "생부모 비동거",
        TRUE ~ NA_character_
      ),
      levels = c("양부모 동거", "모만 동거", "부만 동거", "생부모 비동거")
    )
  )

# 변수별 결측 수를 출력하여 분석 대상에서 제외되는 사례를 점검합니다.
missing_by_var <- dat %>%
  select(all_of(c("book_f", control_vars, pv_lit))) %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "변수", values_to = "결측") %>%
  filter(결측 > 0) %>%
  arrange(desc(결측))

cat("\n[변수별 결측]\n")
print(missing_by_var, n = Inf)

# 모든 PV·분석변수·가중치에 결측이 없는 동일한 표본을 전체 분석에 사용합니다.
# 도서 보유 수는 10권 이하, 부모 교육수준은 중졸 이하를 기준집단으로 고정합니다.
cc_data <- dat %>%
  select(all_of(c(pv_lit, "book_f", control_vars, main_weight, rep_weights))) %>%
  tidyr::drop_na() %>%
  mutate(across(where(is.factor), droplevels)) %>%
  mutate(
    book_f       = relevel(book_f, ref = book_ref),
    parent_edu_f = relevel(parent_edu_f, ref = "중졸 이하"),
    w_           = .data[[main_weight]]
  )

if (any(!is.finite(cc_data[[main_weight]])) || any(cc_data[[main_weight]] <= 0) ||
    any(!is.finite(as.matrix(cc_data[, rep_weights]))) ||
    any(as.matrix(cc_data[, rep_weights]) < 0)) {
  stop("최종가중치는 양수, 반복가중치는 0 이상이며 모두 유한해야 합니다.")
}

# 완전사례 표본으로 설계를 구성하고, PV 통합에 사용할 설계 자유도를 구합니다.
cc_n   <- nrow(cc_data)
cc_des <- make_design(cc_data)
cc_df  <- survey::degf(cc_des)

cat("\n[표본]\n")
cat(sprintf("  원자료           : %s명\n", format(nrow(dat_raw), big.mark = ",")))
cat(sprintf("  결측 제외        : %s명(%.2f%%)\n",
            format(nrow(dat_raw) - cc_n, big.mark = ","),
            100 * (nrow(dat_raw) - cc_n) / nrow(dat_raw)))
cat(sprintf("  최종 분석 표본   : %s명\n", format(cc_n, big.mark = ",")))
cat(sprintf("  복합표본 설계 자유도: %d\n", cc_df))


# 6. 문해력 PV별로 동일한 복합표본 회귀모형을 적합한 뒤 통합합니다.
main_model_list <- map(pv_lit, function(y) {
  survey::svyglm(make_reg_formula(y), design = cc_des)
})
names(main_model_list) <- pv_lit

main_pooled <- pool_svyglm_models(main_model_list, df_com = cc_df, n_used = cc_n)

main_effect_results <- main_pooled %>%
  filter(term %in% SENSE_TREAT_TERMS) %>%
  mutate(`도서 보유 수` = factor(str_remove(term, "^book_f"),
                            levels = setdiff(book_levels, book_ref))) %>%
  arrange(`도서 보유 수`)

# 각 PV의 가중 분산을 평균한 뒤 제곱근을 취해 효과크기(SD)의 분모로 사용합니다.
pv_var_vec <- vapply(pv_lit, function(y) {
  as.numeric(survey::svyvar(as.formula(paste0("~", y)), design = cc_des))
}, numeric(1))

sd_lit <- sqrt(mean(pv_var_vec))

cat(sprintf("\n[문해력 가중 표준편차] %.2f점\n", sd_lit))


# 7. 기술통계: 배경변수별 표기 순서와 범주를 지정합니다.
background_specs <- list(
  list(var = "gender_f",      label = "성별",
       levels = c("남성", "여성")),
  list(var = "age_f",         label = "연령",
       levels = c("16-19세", "20-24세", "25-29세", "30-34세", "35-39세",
                  "40-44세", "45-49세", "50-54세", "55-59세", "60-65세")),
  list(var = "parent_edu_f",  label = "부모 교육수준",
       levels = c("중졸 이하", "고졸 수준", "대졸 이상")),
  list(var = "mother_work_f", label = "모 경제활동",
       levels = c("유급직", "무직/가사 등")),
  list(var = "father_work_f", label = "부 경제활동",
       levels = c("유급직", "무직/가사 등")),
  list(var = "residence_f",   label = "거주지역 규모",
       levels = c("대도시", "중소도시", "소도시/읍면", "농어촌/시골")),
  list(var = "family14_f",    label = "가족구조",
       levels = c("양부모 동거", "모만 동거", "부만 동거", "생부모 비동거"))
)

# 도서 보유 수 범주의 고정 순서에 맞춰 표 1의 한 행을 만듭니다.
make_row <- function(label, values_chr, p_text = "") {
  values <- as.list(values_chr[book_levels])
  names(values) <- unname(book_col_labels[book_levels])
  bind_cols(tibble(변수 = label), as_tibble(values), tibble(`pᵃ` = p_text))
}

# 집단 간 분포 차이는 Rao–Scott F 검정을 먼저 시도합니다.
# 계산에 실패하면 adjWald, Chisq 순으로 대체하며, 모두 실패하면 NA를 반환합니다.
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

# 도서 보유 수 집단별로 최종가중치를 합산하여 범주의 열 백분율을 계산합니다.
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

# 배경변수별 검정 p값과 범주별 가중 백분율을 하나의 표 구획으로 묶습니다.
make_table1_section <- function(spec) {
  empty <- setNames(rep("", length(book_levels)), book_levels)

  header <- make_row(spec$label, empty, fmt_p_with_stars(safe_svy_chisq_p(spec$var)))
  rows   <- map_dfr(spec$levels, function(cat_label) {
    make_row(cat_label, fmt_num(weighted_col_percent(spec$var, cat_label), 1))
  })

  bind_rows(header, rows)
}

# 표 1 첫 두 행: 비가중 사례 수와 최종가중치를 적용한 집단 비율.
book_distribution <- cc_data %>%
  group_by(book_f) %>%
  summarise(n = n(), weighted_n = sum(.data[[main_weight]]), .groups = "drop") %>%
  mutate(pct = 100 * weighted_n / sum(weighted_n))

n_vec   <- setNames(as.character(book_distribution$n),
                    as.character(book_distribution$book_f))
pct_vec <- setNames(fmt_num(book_distribution$pct, 1),
                    as.character(book_distribution$book_f))

표1_기술통계_교차분석 <- bind_rows(
  make_row("n(비가중)", n_vec),
  make_row("%(가중)",   pct_vec),
  map_dfr(background_specs, make_table1_section)
)

show_table(
  표1_기술통계_교차분석,
  "<표 1> 14세 시점 가정 내 도서 보유 수 범주별 주요 배경변수의 가중 분포",
  "table1.csv"
)


# 8. 통합 회귀계수, 표준오차, t값, p값, 95% CI 및 표준화한 차이를 표 2로 출력합니다.
표2_회귀분석 <- main_effect_results %>%
  transmute(
    `도서 보유 수`  = as.character(`도서 보유 수`),
    추정치          = fmt_num(estimate, 2),
    표준오차        = fmt_num(std_error, 2),
    t               = fmt_num(statistic, 2),
    p               = fmt_p_with_stars(p_value),
    `95% 신뢰구간`  = fmt_ci(conf_low, conf_high),
    `효과크기(SD)`  = fmt_num(estimate / sd_lit, 2)
  )

show_table(
  표2_회귀분석,
  "<표 2> 14세 시점 가정 내 도서 보유 수가 성인기 문해력에 미치는 효과 추정 결과",
  "table2.csv"
)


# 9. 민감도 분석
# 지정한 계수 집단을 추가했을 때 감소하는 가중 잔차제곱합으로 집단 부분 R2를 구합니다.
# extra_ss / (extra_ss + full_ss): 축소모형의 잔차 중 해당 집단이 설명하는 비율.
# WLS의 비표준화 공분산을 사용하며, 복합표본 t값을 부분 R2로 변환하지 않습니다.
group_r2_wls <- function(model, terms) {
  sm <- summary(model)
  cf <- coef(model)
  if (!length(terms) || !all(terms %in% names(cf)) || anyNA(cf)) {
    stop("비교기준 계수 이름 또는 모형의 완전계수를 확인하세요.")
  }
  b <- cf[terms]
  v <- sm$cov.unscaled[terms, terms, drop = FALSE]
  extra_ss <- as.numeric(crossprod(b, solve(v, b)))
  full_ss <- sum(model$weights * residuals(model)^2)
  if (!is.finite(extra_ss) || !is.finite(full_ss) || full_ss <= 0 || extra_ss < -1e-8) {
    stop("집단 부분 R2 계산 오류")
  }
  max(0, extra_ss) / (max(0, extra_ss) + full_ss)
}

# 편향 척도와 부분 R2 계산을 위한 WLS 모형 및 범주형 변수의 설계행렬을 준비합니다.
lm_model_list <- setNames(lapply(pv_lit, function(y) {
  stats::lm(make_reg_formula(y), data = cc_data, weights = w_)
}), pv_lit)
mm <- model.matrix(make_reg_formula(pv_lit[1]), data = cc_data)
mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
name_map <- setNames(paste0("v", seq_len(ncol(mm))), colnames(mm))
mm_df <- as.data.frame(mm)
names(mm_df) <- unname(name_map)
mm_df$w_ <- cc_data$w_

# 주 분석의 비교기준은 부모 교육수준의 두 비기준 더미를 묶은 집단입니다.
# 대졸 이상 더미 하나를 사용하는 비교 결과도 함께 출력합니다(원고 표 3은 전체 집단 사용).
benchmark_groups <- list(
  "부모 교육수준 전체" = grep("^parent_edu_f", colnames(mm), value = TRUE),
  "부모 교육수준: 대졸 이상" = SENSE_BENCHMARK_COV
)
if (any(lengths(benchmark_groups) == 0) ||
    !all(unlist(benchmark_groups) %in% colnames(mm))) {
  stop("부모 교육수준의 비교기준 범주가 누락되어 있습니다.")
}

# 통합 계수의 부호를 기준으로 편향 차감 방향을 고정합니다.
# 같은 방향을 모든 PV와 반복표본에 적용하며, 큰 교란에서는 조정치의 부호가 바뀔 수 있습니다.
sense_direction <- setNames(ifelse(main_effect_results$estimate >= 0, 1, -1),
                            main_effect_results$term)

# PV i와 도서 범주 tt에 대해 최종가중치 및 각 반복가중치에서 필요한 양을 저장합니다.
# 각 반복에서 결과모형·원인모형·비교기준 집단의 부분 R2를 모두 다시 계산합니다.
# 원인모형은 해당 도서 범주 더미를 다른 도서 범주 더미와 통제변수에 회귀한 모형입니다.
make_sense_cache <- function(i, tt) {
  message("민감도 반복추정: ", tt, " / ", pv_lit[i])
  s <- sense_direction[[tt]]
  survey::withReplicates(
    design = cc_des,
    theta = function(w, data) {
      if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
        stop("분석할 수 없는 반복가중치입니다.")
      }
      data$w_ <- as.numeric(w)
      m <- stats::lm(make_reg_formula(pv_lit[i]), data = data, weights = w_)
      tab <- coef(summary(m))
      if (anyNA(coef(m)) || !tt %in% rownames(tab)) stop("결과모형 특이행렬: ", tt)
      beta <- unname(tab[tt, "Estimate"])
      # a는 WLS 계수의 SE × sqrt(잔차 자유도)로, 편향식에 쓰는 잔차 척도 비율입니다.
      # 최종 추론의 표준오차는 아래 반복분산과 PV 통합으로 따로 계산합니다.
      a <- unname(tab[tt, "Std. Error"] * sqrt(m$df.residual))
      if (!is.finite(a) || a <= 0) stop("편향 척도가 양수가 아닙니다.")

      ddata <- mm_df
      ddata$w_ <- as.numeric(w)
      y_nm <- name_map[[tt]]
      x_nm <- setdiff(unname(name_map), y_nm)
      md <- stats::lm(stats::reformulate(x_nm, response = y_nm),
                      data = ddata, weights = w_)
      if (anyNA(coef(md))) stop("노출모형 특이행렬: ", tt)
      out <- c(beta = beta, bias_scale = a, wls_df = m$df.residual,
               wls_t = unname(tab[tt, "t value"]))
      for (g in seq_along(benchmark_groups)) {
        terms <- benchmark_groups[[g]]
        ry <- group_r2_wls(m, terms)
        rd <- group_r2_wls(md, unname(name_map[terms]))
        out[paste0("rd_", g)] <- rd
        out[paste0("ry_", g)] <- ry
        for (j in seq_along(SENSE_KD)) {
          # kd=ky에 1·2·3을 적용하여 비교기준에 상대적인 교란 강도 상한을 구합니다.
          bounds <- sensemakr::ovb_partial_r2_bound(
            r2dxj.x = rd, r2yxj.dx = ry,
            kd = SENSE_KD[j], ky = SENSE_KD[j]
          )
          r_d <- bounds$r2dz.x
          r_y <- bounds$r2yz.dx
          if (length(r_d) != 1L || length(r_y) != 1L ||
              any(!is.finite(c(r_d, r_y))) ||
              r_d < 0 || r_d >= 1 || r_y < 0 || r_y > 1) {
            stop("비교기준 시나리오가 허용 범위를 벗어났습니다: ", tt,
                 " / ", names(benchmark_groups)[g], " / ", SENSE_KD[j], "배")
          }
          nm <- paste0(g, "_", j)
          # 조정치 = beta - s × a × sqrt(rY × rD / (1-rD)).
          # rD: 나머지 설명변수 통제 후 미관측 교란의 원인변수 부분설명력.
          # rY: 원인변수와 나머지 설명변수 통제 후 미관측 교란의 결과변수 부분설명력.
          out[paste0("adj_", nm)] <- beta - s*a*sqrt(r_y*r_d/(1-r_d))
          out[paste0("bound_d_", nm)] <- r_d
          out[paste0("bound_y_", nm)] <- r_y
        }
      }
      if (any(!is.finite(out))) stop("민감도 계산에 비유한 값이 있습니다.")
      out
    },
    return.replicates = TRUE
  )
}

sense_cache <- setNames(lapply(SENSE_TREAT_TERMS, function(tt) {
  lapply(seq_along(pv_lit), function(i) make_sense_cache(i, tt))
}), SENSE_TREAT_TERMS)

# 저장된 최종가중치 통계량과 반복추정치를 재사용하여 임계값 탐색의 중복 적합을 줄입니다.
cache_stats <- setNames(lapply(SENSE_TREAT_TERMS, function(tt) {
  lapply(sense_cache[[tt]], function(obj) {
    z <- coef(obj)
    reps <- as.matrix(obj$replicates)
    if (ncol(reps) != length(z)) stop("반복추정치 열 수 불일치")
    colnames(reps) <- names(z)
    list(z = z, reps = reps)
  })
}), SENSE_TREAT_TERMS)

# PV별 조정치와 반복분산을 반환합니다. f = sqrt(rY × rD / (1-rD)).
# j=NULL: RV 탐색용으로 같은 f(따라서 공통 rD, rY 경로)를 PV와 반복표본에 적용합니다.
# j 지정: g번째 비교기준·j번째 배수에서 반복표본별로 재계산한 조정치를 사용합니다.
# 반복분산은 최종가중치 조정치를 중심으로 계산하며, 단순 SE 배율 근사가 아닙니다.
pv_sense_rows <- function(tt, f = 0, j = NULL, g = 1L) {
  s <- sense_direction[[tt]]
  bind_rows(lapply(seq_along(pv_lit), function(i) {
    x <- cache_stats[[tt]][[i]]
    if (is.null(j)) {
      q <- unname(x$z["beta"] - s*f*x$z["bias_scale"])
      qr <- x$reps[, "beta"] - s*f*x$reps[, "bias_scale"]
    } else {
      nm <- paste0("adj_", g, "_", j)
      q <- unname(x$z[nm])
      qr <- x$reps[, nm]
    }
    u <- as.numeric(survey::svrVar(
      thetas = qr, scale = cc_des$scale, rscales = cc_des$rscales,
      mse = TRUE, coef = q
    ))
    if (!is.finite(q) || !is.finite(u) || u < 0) stop("유효하지 않은 조정량/반복분산")
    tibble(treatment = tt, pv = pv_lit[i], est = q, se = sqrt(u))
  }))
}

# PV별 조정치와 반복분산을 결합하고, SENSE_ALPHA에 따른 CI와 양측 p값을 계산합니다.
pool_sense <- function(tt, f = 0, j = NULL, g = 1L) {
  ans <- pool_scalar_df(
    df = pv_sense_rows(tt, f, j, g), key_cols = "treatment",
    estimate_col = "est", se_col = "se", df_com = cc_df, n_used = cc_n
  )
  crit <- qt(1-SENSE_ALPHA/2, df = ans$df)
  ans$conf_low <- ans$estimate - crit*ans$std_error
  ans$conf_high <- ans$estimate + crit*ans$std_error
  ans$p_value <- 2*pt(-abs(ans$statistic), df = ans$df)
  ans$stars <- stars_from_p(ans$p_value)
  ans
}

# 교란 강도 0에서 주 회귀분석의 계수·SE·자유도·CI·p값이 재현되는지 확인합니다.
for (tt in SENSE_TREAT_TERMS) {
  z <- pool_sense(tt)
  ref <- main_effect_results[main_effect_results$term == tt, ]
  for (nm in c("estimate", "std_error", "df", "conf_low", "conf_high", "p_value")) {
    if (!isTRUE(all.equal(z[[nm]], ref[[nm]], tolerance = 1e-6))) {
      stop("교란 0 원모형 재현 실패: ", tt, " / ", nm)
    }
  }
}
message("교란 0에서 계수·SE·자유도·CI·p값 재현: 통과")

# 모든 도서 범주 × 비교기준 × 배수 조합의 통합 조정 결과를 구성합니다.
sensitivity_bounds_pooled <- bind_rows(lapply(SENSE_TREAT_TERMS, function(tt) {
  bind_rows(lapply(seq_along(benchmark_groups), function(g) {
    bind_rows(lapply(seq_along(SENSE_KD), function(j) {
      z <- pool_sense(tt, j = j, g = g)
      z$benchmark <- names(benchmark_groups)[g]
      z$scenario <- scenario_levels[j]
      z
    }))
  }))
}))

# 조정된 통합 CI가 처음 0을 포함하는 편향계수 f의 경계를 수치적으로 구합니다.
# 일반 회귀용 sensemakr의 유의성 기준 RV 닫힌식이 아닌, 이 파일의 통합 CI를 사용합니다.
# margin은 원래 계수의 부호 방향에서 0에 가까운 CI 끝점이며, 0 이하가 되면 유의성이 사라집니다.
find_sense_threshold <- function(tt) {
  s <- sense_direction[[tt]]
  margin <- function(f) {
    z <- pool_sense(tt, f = f)
    s*z$estimate - qt(1-SENSE_ALPHA/2, z$df)*z$std_error
  }
  if (margin(0) <= 0) return(0)
  z0 <- pool_sense(tt)
  a_bar <- mean(vapply(cache_stats[[tt]], function(x) unname(x$z["bias_scale"]), numeric(1)))
  f_zero <- abs(z0$estimate)/a_bar
  # 0부터 통합 점추정치가 0이 되는 f까지 격자를 탐색한 뒤 첫 교차구간에서 근을 구합니다.
  # 두 격자 해상도(129, 257)에서 얻은 임계값을 비교하여 수치 안정성을 확인합니다.
  solve_grid <- function(grid_points) {
    grid <- seq(0, f_zero, length.out = grid_points)
    vals <- vapply(grid, margin, numeric(1))
    hit <- which(vals <= 0)[1]
    if (is.na(hit)) {
      # 점추정치 영점에서 발생할 수 있는 작은 부동소수점 오차만 허용합니다.
      if (abs(tail(vals, 1)) < 1e-9) return(f_zero)
      stop("CI 임계값을 찾지 못했습니다: ", tt)
    }
    if (hit == 1L) return(0)
    uniroot(margin, grid[c(hit-1L, hit)], tol = 1e-10)$root
  }
  f1 <- solve_grid(129L)
  f2 <- solve_grid(257L)
  if (abs(f1-f2) > 1e-6*max(1, f2)) stop("임계값 탐색 해상도 불안정: ", tt)
  f2
}

# f*를 부분설명력 척도로 변환합니다.
# RV05: rD=rY=r일 때 f*=r/sqrt(1-r)를 만족하는 r.
# extreme_RV05: 0<=rY<=1을 허용할 때 필요한 최소 rD = f*^2/(1+f*^2).
# rY=1은 잔여 결과변량을 전부 설명하는 형식적 경계입니다.
# 임계점의 p값, 직전 p값 및 CI를 함께 저장하여 탐색 결과를 점검합니다.
sensitivity_thresholds <- bind_rows(lapply(SENSE_TREAT_TERMS, function(tt) {
  message("유의성 소실 임계값 탐색: ", tt)
  f_star <- find_sense_threshold(tt)
  rv <- if (f_star == 0) 0 else 2/(1+sqrt(1+4/f_star^2))
  xrv <- f_star^2/(1+f_star^2)
  z <- pool_sense(tt, f_star)
  zb <- pool_sense(tt, max(0, f_star-1e-5*max(1, f_star)))
  if (f_star > 0 && abs(z$p_value-SENSE_ALPHA) > 1e-5) {
    stop("임계점에서 p값이 alpha와 일치하지 않습니다: ", tt)
  }
  # 부분설명력 변환의 역산과 extreme_RV05 <= RV05 관계를 확인합니다.
  if (f_star > 0 && (abs(rv/sqrt(1-rv)-f_star)>1e-7 ||
                    abs(sqrt(xrv/(1-xrv))-f_star)>1e-7)) stop("RV 변환 점검 실패")
  if (xrv > rv+1e-9) stop("XRV가 RV보다 큽니다.")
  tibble(treatment = tt, RV05 = rv, extreme_RV05 = xrv,
         f_threshold = f_star, estimate_at_threshold = z$estimate,
         se_at_threshold = z$std_error, df_at_threshold = z$df,
         p_at_threshold = z$p_value, p_just_below = zb$p_value,
         conf_low_at_threshold = z$conf_low,
         conf_high_at_threshold = z$conf_high)
}))

# PV별 WLS t값·자유도, 편향 척도, 비교기준 부분 R2 및 교란 상한을 재현용으로 저장합니다.
sensitivity_pv_audit <- bind_rows(lapply(SENSE_TREAT_TERMS, function(tt) {
  bind_rows(lapply(seq_along(pv_lit), function(i) {
    z <- cache_stats[[tt]][[i]]$z
    bind_rows(lapply(seq_along(benchmark_groups), function(g) {
      bind_rows(lapply(seq_along(SENSE_KD), function(j) {
        nm <- paste0(g,"_",j)
        tibble(treatment = tt, pv = pv_lit[i],
               benchmark = names(benchmark_groups)[g], scenario = scenario_levels[j],
               beta = unname(z["beta"]), bias_scale = unname(z["bias_scale"]),
               wls_t = unname(z["wls_t"]), wls_df = unname(z["wls_df"]),
               benchmark_rD = unname(z[paste0("rd_",g)]),
               benchmark_rY = unname(z[paste0("ry_",g)]),
               bound_rD = unname(z[paste0("bound_d_",nm)]),
               bound_rY = unname(z[paste0("bound_y_",nm)]))
      }))
    }))
  }))
}))

# 표 3: 두 비교기준의 조정치(SE) [95% CI]와 RV를 함께 출력합니다.
# RV는 비교기준과 무관한 임계값이므로 두 비교기준 행에서 동일하게 표시됩니다.
표3_민감도분석 <- sensitivity_bounds_pooled %>%
  mutate(benchmark = factor(benchmark, levels = names(benchmark_groups)),
         scenario = factor(scenario, levels = scenario_levels),
         value = paste0(fmt_num(estimate,2), stars_from_p(p_value),
                        " (",fmt_num(std_error,2),") ",fmt_ci(conf_low, conf_high))) %>%
  select(treatment, benchmark, scenario, value) %>%
  pivot_wider(names_from = scenario, values_from = value) %>%
  left_join(sensitivity_thresholds %>% select(treatment, RV05, extreme_RV05), by="treatment") %>%
  arrange(match(treatment, SENSE_TREAT_TERMS), benchmark) %>%
  transmute(`도서 보유 수` = str_remove(treatment,"^book_f"),
            `비교기준` = as.character(benchmark),
            RV05 = fmt_dec(RV05,3), `극단적 RV05` = fmt_dec(extreme_RV05,3),
            across(all_of(scenario_levels)))
show_table(표3_민감도분석, "<표 3> RV05·극단적 RV05 및 조정치(SE) [95% CI]", "table3.csv")
show_table(sensitivity_bounds_pooled, "[상세] 조정 추정치·SE·df·CI·p값", "sensitivity_details.csv")
show_table(sensitivity_thresholds, "[점검] 유의성 소실 임계값과 경계 p값", "threshold_audit.csv")
show_table(sensitivity_pv_audit, "[재현] PV별 편향 척도 및 비교기준 부분 R2", "sensitivity_pv_audit.csv")

cat("\n[해석] RV05와 극단적 RV05는 통합 95% CI가 최초로 0을 포함하는 교란 강도입니다.\n")
cat("[주의] 이 값들은 일반 OLS용 sensemakr alpha=.05의 닫힌식을 그대로 호출한 값이 아닙니다.\n")
cat("[자유도] 완전자료 자유도는 survey::degf(cc_des), PV 통합은 Barnard-Rubin 형태의 소표본 보정입니다.\n")
cat("[주의] rY=1은 극단적인 형식적 경계이며, 실제로 그런 교란이 존재한다는 의미가 아닙니다.\n")
cat("[주의] 비교기준 시나리오는 반복표본별 재보정, RV는 고정 부분 R2 경로를 사용합니다.\n")

# 상세 객체와 sessionInfo를 저장하여 추정치 및 실행 환경을 확인할 수 있게 합니다.
if (isTRUE(SAVE_OUTPUT)) {
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(main_results=main_effect_results,
               adjusted_results=sensitivity_bounds_pooled,
               thresholds=sensitivity_thresholds, pv_audit=sensitivity_pv_audit,
               cache=sense_cache, alpha=SENSE_ALPHA, benchmarks=benchmark_groups,
               design_df=cc_df, fay_rho=fay_rho, session=sessionInfo()),
          file.path(OUTPUT_DIR,"sensitivity_results.rds"))
  writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR,"session_info.txt"))
}

# 10. 실행 환경: R 및 주요 패키지 버전을 콘솔에 출력합니다.
cat("\n[재현 정보]\n", R.version.string, "\n", sep="")
for (pkg in packages) cat(sprintf("   %-12s %s\n", pkg, as.character(packageVersion(pkg))))
