# ==============================================================================
# 청소년기 가정 내 도서 보유 수가 성인기 문해력에 미치는 영향
# PIAAC 2주기 한국 자료 분석 코드
#
# Data : PIAAC Cycle 2 Korea public use file (prgkorp2.csv, OECD)
#        자료는 저장소에 포함하지 않았다. OECD에서 내려받아 이 스크립트와
#        같은 폴더에 두면 그대로 실행된다.
#        배포본에 따라 구분자가 세미콜론(;), 쉼표(,), 탭인 경우가 있어
#        첫 줄에서 자동으로 판별한다.
#
# 산출물
#   표 1. 도서 보유 수 범주별 배경변수의 가중 분포 및 교차분석
#   표 2. 도서 보유 수와 성인기 문해력의 회귀분석(효과크기 포함)
#   표 3. 미관측 교란요인에 대한 민감도 분석
#
# 분석 원칙
#   - 결과변수는 문해력 유의측정값(PVLIT1-PVLIT10)이다.
#   - 최종가중치(SPFWT0)와 80개 반복가중치(SPFWT1-SPFWT80)를 적용한다.
#   - 분산 추정은 Fay의 균형반복복제법(BRR)이며 rho = 0.3이다.
#   - PV별 복합표본 추정치는 Rubin의 결합규칙으로 통합한다.
#
# 구성
#   0 실행 옵션      1 패키지        2 분석 설정      3 유틸리티 함수
#   4 결합규칙 함수  5 자료 전처리   6 모형 적합      7 표 1
#   8 표 2           9 표 3          10 재현 정보
# ==============================================================================


# ==============================================================================
# 0. 실행 옵션
# ==============================================================================

DATA_PATH   <- "prgkorp2.csv"
USE_VIEWER  <- interactive()   # RStudio 표 창 자동 실행 여부
SAVE_OUTPUT <- FALSE           # TRUE이면 표를 CSV로 저장
OUTPUT_DIR  <- "output"


# ==============================================================================
# 1. 패키지
# ==============================================================================

packages <- c("tidyverse", "survey", "sensemakr")

to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

invisible(lapply(packages, library, character.only = TRUE))


# ==============================================================================
# 2. 분석 설정
# ==============================================================================

## (1) 가중치 및 결과변수 -----------------------------------------------------

pv_lit      <- paste0("PVLIT", 1:10)   # 문해력 유의측정값
main_weight <- "SPFWT0"                # 최종가중치
rep_weights <- paste0("SPFWT", 1:80)   # 반복가중치
fay_rho     <- 0.3                     # Fay 조정계수

## (2) 원인변수 및 통제변수 ---------------------------------------------------

book_var <- "J2_Q06"      # 14세 시점 가정 내 도서 보유 수
age_var  <- "AGEG5LFS"    # 연령 범주

book_levels <- c(
  "10권 이하", "11-25권", "26-100권",
  "101-200권", "201-500권", "500권 초과"
)
book_ref <- "10권 이하"   # 기준집단

book_col_labels <- c(
  "10권 이하"  = "≤10",
  "11-25권"    = "11-25",
  "26-100권"   = "26-100",
  "101-200권"  = "101-200",
  "201-500권"  = "201-500",
  "500권 초과" = ">500"
)

control_vars <- c(
  "gender_f",       # 성별
  "age_f",          # 연령
  "parent_edu_f",   # 부모 교육수준
  "mother_work_f",  # 14세 당시 모 경제활동
  "father_work_f",  # 14세 당시 부 경제활동
  "residence_f",    # 14세 당시 거주지역 규모
  "family14_f"      # 14세 당시 가족구조
)

## (3) 민감도 분석 설정 -------------------------------------------------------

# 기준집단(10권 이하)과 비교되는 모든 도서 보유 수 범주
SENSE_TREAT_TERMS   <- paste0("book_f", setdiff(book_levels, book_ref))
SENSE_BENCHMARK_COV <- "parent_edu_f대졸 이상"     # 비교기준 공변량
SENSE_KD            <- c(1, 2, 3)                # 비교기준 공변량 대비 배수
SENSE_ALPHA         <- 0.05                      # 강건성 판정 유의수준

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
    rho              = fay_rho,
    combined.weights = TRUE
  )
}

# 회귀식: 문해력 PV ~ 도서 보유 수 + 통제변수
make_reg_formula <- function(outcome, controls = control_vars) {
  as.formula(paste(outcome, "~ book_f +", paste(controls, collapse = " + ")))
}

## 출력 형식 ------------------------------------------------------------------

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

# 앞자리 0을 생략한 소수 표기(p, R² 등)
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

# 표 주석 출력(줄바꿈 정리)
show_note <- function(...) {
  cat("\n", paste(strwrap(paste0(...), width = 84), collapse = "\n"), "\n", sep = "")
}


# ==============================================================================
# 4. Rubin 결합규칙
# ==============================================================================

## (1) svyglm 모형 목록(PV별)의 회귀계수 통합 ---------------------------------

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

## (2) 단일 추정치·표준오차 형태(민감도 분석 결과 등)의 통합 -------------------
##     key_cols는 하나 이상 지정해야 한다.

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


# ==============================================================================
# 5. 자료 불러오기 및 전처리
# ==============================================================================

## (1) 구분자 자동 판별 후 읽기 -----------------------------------------------

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

## (2) 필요한 변수 확인 -------------------------------------------------------

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

## (3) 변수 생성 --------------------------------------------------------------

special_numeric_vars <- c(
  book_var, "PAREDC2", "J2_Q04d", "J2_Q05d",
  "J2_Q07_C", "J2_Q0801", "J2_Q0802"
)

dat <- dat_raw %>%
  mutate(across(where(is.character), recode_char_missing)) %>%
  mutate(across(all_of(special_numeric_vars), recode_numeric_special_missing)) %>%
  mutate(across(all_of(c(pv_lit, main_weight, rep_weights)),
                ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(
    # 14세 시점 가정 내 도서 보유 수(기준범주가 첫 수준)
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

    # 성별
    gender_f = factor(
      case_when(
        as_code(GENDER_R) == 1 ~ "남성",
        as_code(GENDER_R) == 2 ~ "여성",
        TRUE ~ NA_character_
      ),
      levels = c("남성", "여성")
    ),

    # 연령
    age_f = factor(
      case_when(
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
      levels = c("20-24세", "25-29세", "30-34세", "35-39세", "40-44세",
                 "45-49세", "50-54세", "55-59세", "60-64세", "65세 이상")
    ),

    # 부모 교육수준(기준범주가 첫 수준)
    parent_edu_f = factor(
      case_when(
        as_code(PAREDC2) == 1 ~ "중졸 이하",
        as_code(PAREDC2) == 2 ~ "고졸/전문대졸",
        as_code(PAREDC2) == 3 ~ "대졸 이상",
        TRUE ~ NA_character_
      ),
      levels = c("중졸 이하", "고졸/전문대졸", "대졸 이상")
    ),

    # 14세 당시 모 경제활동
    mother_work_f = factor(
      case_when(
        as_code(J2_Q04d) == 1 ~ "유급직",
        as_code(J2_Q04d) == 2 ~ "무직/가사 등",
        TRUE ~ NA_character_
      ),
      levels = c("유급직", "무직/가사 등")
    ),

    # 14세 당시 부 경제활동
    father_work_f = factor(
      case_when(
        as_code(J2_Q05d) == 1 ~ "유급직",
        as_code(J2_Q05d) == 2 ~ "무직/가사 등",
        TRUE ~ NA_character_
      ),
      levels = c("유급직", "무직/가사 등")
    ),

    # 14세 당시 거주지역 규모
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

    # 14세 당시 가족구조
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

## (4) 변수별 결측 점검 -------------------------------------------------------

missing_by_var <- dat %>%
  select(all_of(c("book_f", control_vars, pv_lit))) %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "변수", values_to = "결측") %>%
  filter(결측 > 0) %>%
  arrange(desc(결측))

cat("\n[변수별 결측]\n")
print(missing_by_var, n = Inf)

## (5) 완전사례 자료 ----------------------------------------------------------

cc_data <- dat %>%
  select(all_of(c(pv_lit, "book_f", control_vars, main_weight, rep_weights))) %>%
  tidyr::drop_na() %>%
  mutate(across(where(is.factor), droplevels)) %>%
  mutate(
    book_f       = relevel(book_f, ref = book_ref),
    parent_edu_f = relevel(parent_edu_f, ref = "중졸 이하"),
    w_           = .data[[main_weight]]          # lm의 weights 인자용
  )

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


# ==============================================================================
# 6. 모형 적합 및 결합
# ==============================================================================

## (1) PV별 복합표본 회귀 ------------------------------------------------------

main_model_list <- map(pv_lit, function(y) {
  survey::svyglm(make_reg_formula(y), design = cc_des)
})
names(main_model_list) <- pv_lit

## (2) Rubin 결합 및 원인변수 계수 추출 ---------------------------------------

main_pooled <- pool_svyglm_models(main_model_list, df_com = cc_df, n_used = cc_n)

main_effect_results <- main_pooled %>%
  filter(term %in% SENSE_TREAT_TERMS) %>%
  mutate(`도서 보유 수` = factor(str_remove(term, "^book_f"),
                            levels = setdiff(book_levels, book_ref))) %>%
  arrange(`도서 보유 수`)

## (3) 문해력의 가중 표준편차(효과크기 환산 기준) -----------------------------
##     PV별 가중 분산의 평균(Rubin 결합의 점추정치)에 제곱근을 취한다.

pv_var_vec <- vapply(pv_lit, function(y) {
  as.numeric(survey::svyvar(as.formula(paste0("~", y)), design = cc_des))
}, numeric(1))

sd_lit <- sqrt(mean(pv_var_vec))

cat(sprintf("\n[문해력 가중 표준편차] %.2f점\n", sd_lit))


# ==============================================================================
# 7. 표 1. 도서 보유 수 범주별 배경변수의 가중 분포 및 교차분석
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
  list(var = "residence_f",   label = "거주지역 규모",
       levels = c("대도시", "중소도시", "소도시/읍면", "농어촌/시골")),
  list(var = "family14_f",    label = "가족구조",
       levels = c("양부모 동거", "모만 동거", "부만 동거", "생부모 비동거"))
)

# 표 1의 한 행 만들기(값은 이미 문자형으로 정리된 벡터)
make_row <- function(label, values_chr, p_text = "") {
  values <- as.list(values_chr[book_levels])
  names(values) <- unname(book_col_labels[book_levels])
  bind_cols(tibble(변수 = label), as_tibble(values), tibble(`pᵃ` = p_text))
}

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

make_table1_section <- function(spec) {
  empty <- setNames(rep("", length(book_levels)), book_levels)

  header <- make_row(spec$label, empty, fmt_p_with_stars(safe_svy_chisq_p(spec$var)))
  rows   <- map_dfr(spec$levels, function(cat_label) {
    make_row(cat_label, fmt_num(weighted_col_percent(spec$var, cat_label), 1))
  })

  bind_rows(header, rows)
}

## (1) 범주별 사례수 및 가중 구성비 -------------------------------------------

book_distribution <- cc_data %>%
  group_by(book_f) %>%
  summarise(n = n(), weighted_n = sum(.data[[main_weight]]), .groups = "drop") %>%
  mutate(pct = 100 * weighted_n / sum(weighted_n))

n_vec   <- setNames(as.character(book_distribution$n),
                    as.character(book_distribution$book_f))
pct_vec <- setNames(fmt_num(book_distribution$pct, 1),
                    as.character(book_distribution$book_f))

## (2) 표 조립 ----------------------------------------------------------------

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

show_note(
  "주. n(비가중)은 각 범주의 사례수, %(가중)은 최종가중치를 적용한 범주별 ",
  "구성비이다. 배경변수의 수치는 복합표본 최종가중치를 적용한 열 백분율(%)이다. ",
  "a Rao-Scott 수정 카이제곱 검정의 p값이다. ***p < .001."
)


# ==============================================================================
# 8. 표 2. 도서 보유 수와 성인기 문해력의 회귀분석
# ==============================================================================

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

show_note(
  sprintf("주. 표본 크기는 %s명이며, 기준집단은 '%s'이다. ",
          format(cc_n, big.mark = ","), book_ref),
  "모형은 성별, 연령, 부모 교육수준, 14세 당시 모·부의 경제활동 여부, ",
  "거주지역 규모 및 가족구조를 통제하였다. 각 계수와 표준오차는 10개의 문해력 PV에 ",
  sprintf("대해 최종가중치와 80개의 반복가중치, Fay 조정계수 %.1f을 적용한 BRR로 추정한 뒤 ",
          fay_rho),
  "Rubin의 결합규칙에 따라 통합하였다. ",
  sprintf("효과크기(SD)는 각 추정치를 본 표본의 문해력 가중 표준편차(%.2f점)로 나눈 값이다. ",
          sd_lit),
  "***p < .001."
)


# ==============================================================================
# 9. 표 3. 미관측 교란요인에 대한 민감도 분석
# ==============================================================================

# sensemakr는 lm 객체를 요구하므로 부분설명력과 편의 환산척도는
# 최종가중치를 적용한 가중 선형회귀에서 산출한다. 통계적 불확실성은 Fay BRR과
# Rubin의 결합규칙을 적용한 표준오차와 자유도를 통해 반영한다.

lm_model_list <- map(pv_lit, function(y) {
  stats::lm(make_reg_formula(y), data = cc_data, weights = w_)
})
names(lm_model_list) <- pv_lit

## (1) 비교기준 공변량의 관찰된 부분설명력 ------------------------------------
##     표 3 주석에 그대로 들어가는 값이므로 표보다 먼저 계산한다.

# ① 결과변수 쪽: R²_{Y~Z | D, X}
r2_y_benchmark <- mean(vapply(lm_model_list, function(m) {
  as.numeric(sensemakr::partial_r2(m, covariates = SENSE_BENCHMARK_COV))
}, numeric(1)))

# ② 원인변수 쪽: R²_{D~Z | X}
#    sensemakr 내부와 동일하게, 해당 처치 더미를 결과로 두고 나머지 모든 회귀항
#    (다른 도서 보유 수 더미 포함)을 설명변수로 하는 회귀를 적합한다.
#    한글 변수명이 수식에서 문제가 되므로 임시 이름으로 치환한다.
mm <- model.matrix(make_reg_formula(pv_lit[1]), data = cc_data)
mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]

name_map     <- setNames(paste0("v", seq_len(ncol(mm))), colnames(mm))
mm_df        <- as.data.frame(mm)
names(mm_df) <- unname(name_map)
mm_df$w_     <- cc_data$w_

r2_d_benchmark <- vapply(SENSE_TREAT_TERMS, function(tt) {
  y_nm <- name_map[[tt]]
  x_nm <- setdiff(unname(name_map), y_nm)
  m_d  <- stats::lm(as.formula(paste(y_nm, "~", paste(x_nm, collapse = " + "))),
                    data = mm_df, weights = w_)
  as.numeric(sensemakr::partial_r2(m_d, covariates = name_map[[SENSE_BENCHMARK_COV]]))
}, numeric(1))

기준공변량_부분설명력 <- tibble(
  `도서 보유 수`  = str_remove(SENSE_TREAT_TERMS, "^book_f"),
  `R2(원인변수)`  = fmt_dec(r2_d_benchmark, 3),
  `R2(결과변수)`  = fmt_dec(r2_y_benchmark, 3)
)

show_table(
  기준공변량_부분설명력,
  "[보조] 비교기준 공변량('부모 교육수준: 대졸 이상')의 관찰된 부분설명력",
  "benchmark_r2.csv"
)

## (2) 비교기준 공변량 대비 배수 시나리오(bounds) -----------------------------
##     편의 환산척도 se*sqrt(dof)는 가중 lm의 se와 잔차 자유도를 사용하고,
##     조정 표준오차는 PV별 BRR 표준오차에 sensemakr의 조정 비율을 적용한다.
##     이후 PV별 조정 추정치와 분산을 Rubin의 결합규칙으로 통합한다.

benchmark_bounds_pv <- function(i, treatment_term) {
  m       <- lm_model_list[[i]]
  tab_lm  <- coef(summary(m))
  tab_brr <- coef(summary(main_model_list[[i]]))

  if (!(treatment_term %in% rownames(tab_lm))) {
    stop("벤치마크 분석 대상 계수가 모형에 없습니다: ", treatment_term)
  }

  est_lm  <- tab_lm[treatment_term, "Estimate"]
  se_lm   <- tab_lm[treatment_term, "Std. Error"]
  dof_lm  <- m$df.residual
  se_brr  <- tab_brr[treatment_term, "Std. Error"]

  sense_out <- sensemakr::sensemakr(
    model                = m,
    treatment            = treatment_term,
    benchmark_covariates = SENSE_BENCHMARK_COV,
    kd                   = SENSE_KD,
    ky                   = SENSE_KD,
    q                    = 1,
    alpha                = SENSE_ALPHA
  )

  bounds <- as_tibble(sense_out$bounds)
  stopifnot(nrow(bounds) == length(SENSE_KD))

  bounds %>%
    mutate(
      treatment = treatment_term,
      pv        = pv_lit[i],
      scenario  = scenario_levels,
      Adjusted_Estimate = sensemakr::adjusted_estimate(
        estimate = est_lm, se = se_lm, dof = dof_lm,
        r2dz.x = r2dz.x, r2yz.dx = r2yz.dx
      ),
      Adjusted_SE = sensemakr::adjusted_se(
        se = se_brr,
        dof = dof_lm,
        r2dz.x = r2dz.x,
        r2yz.dx = r2yz.dx
      )
    ) %>%
    select(treatment, pv, scenario, r2dz.x, r2yz.dx, Adjusted_Estimate, Adjusted_SE)
}

sensitivity_bounds_all <- map_dfr(SENSE_TREAT_TERMS, function(tt) {
  map_dfr(seq_along(pv_lit), benchmark_bounds_pv, treatment_term = tt)
})

sensitivity_bounds_pooled <- pool_scalar_df(
  df           = sensitivity_bounds_all,
  key_cols     = c("treatment", "scenario"),
  estimate_col = "Adjusted_Estimate",
  se_col       = "Adjusted_SE",
  df_com       = cc_df,
  n_used       = cc_n
)


## (3) 설계·PV 결합을 반영한 RV와 극단적 RV ----------------------------------
##     sensemakr의 부분 R² 모수화에 필요한 회귀 기하량과 BRR·Rubin 결합
##     추론에 필요한 표준오차·자유도를 분리한다. bias_scale은
##     SE_WLS*sqrt(df_WLS)로, 미관측 교란의 부분 R²를 회귀계수 편의로 환산한다.

rv_geometry_pv <- map_dfr(SENSE_TREAT_TERMS, function(tt) {

  map_dfr(seq_along(pv_lit), function(i) {

    m   <- lm_model_list[[i]]
    tab <- coef(summary(m))

    tibble(
      treatment  = tt,
      pv         = pv_lit[i],
      bias_scale = tab[tt, "Std. Error"] * sqrt(m$df.residual),
      dof_geom   = m$df.residual
    )
  })
})

## 공통 부분 R² 시나리오에서 PV별 조정계수를 평균내면
## 편향 환산척도 역시 PV 간 평균으로 결합된다.

rv_geometry_pooled <- rv_geometry_pv %>%
  group_by(treatment) %>%
  summarise(
    bias_scale = mean(bias_scale),
    dof_geom   = mean(dof_geom),
    .groups    = "drop"
  )

## BRR 표본내분산과 PV 간 분산이 이미 결합된 표 2 결과에
## 회귀 기하량을 결합한다.

sensitivity_rv_pooled <- main_effect_results %>%
  transmute(
    treatment       = term,
    `도서 보유 수` = as.character(`도서 보유 수`),
    estimate         = estimate,
    std_error        = std_error,
    t_pooled         = statistic,
    df_pooled        = df
  ) %>%
  left_join(rv_geometry_pooled, by = "treatment") %>%
  mutate(
    ## 미관측 교란 하나를 추가하면 잔차 자유도가 1 감소한다.
    df_test = pmax(df_pooled - 1, 1),

    ## BRR·Rubin 결합 추론의 임계값
    t_critical = qt(
      1 - SENSE_ALPHA / 2,
      df = df_test
    ),

    ## q = 1: 추정효과를 100% 감소시키는 편향 크기
    fq = abs(estimate) / bias_scale,

    ## 통계적 유의성을 제거하는 데 필요한 임계 편향
    ##
    ## sqrt(dof_geom / (dof_geom - 1))은 미관측 변수 하나가
    ## 모형에 추가될 때의 잔차 자유도 변화를 반영한다.
    f_critical =
      t_critical *
      std_error *
      sqrt(dof_geom / (dof_geom - 1)) /
      bias_scale,

    f_difference = pmax(fq - f_critical, 0),

    ## 원인변수·결과변수 부분 R²를 동일하게 둔 RV(α)
    RV = if_else(
      f_difference <= 0,
      0,
      2 / (
        1 + sqrt(
          1 + 4 / (f_difference^2)
        )
      )
    ),

    ## 결과변수 쪽 관련성을 극단적으로 허용한 XRV(α)
    XRV = pmax(
      (fq^2 - f_critical^2) / (1 + fq^2),
      0
    )
  ) %>%
  select(
    treatment,
    `도서 보유 수`,
    estimate,
    std_error,
    t_pooled,
    df_pooled,
    bias_scale,
    RV,
    XRV
  )

cat("\n[설계·PV 결합 강건성 값]\n")

print(
  sensitivity_rv_pooled %>%
    transmute(
      `도서 보유 수`,
      추정치          = round(estimate, 3),
      표준오차        = round(std_error, 3),
      t               = round(t_pooled, 3),
      `결합 자유도`   = round(df_pooled, 1),
      `편향 환산척도` = round(bias_scale, 3),
      `RV(α = .05)`   = round(RV, 3),
      `XRV(α = .05)`  = round(XRV, 3)
    ),
  n = Inf
)

## (4) 표 조립 ----------------------------------------------------------------

표3_벤치마크 <- sensitivity_bounds_pooled %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    조정치   = paste0(fmt_num(estimate, 2), stars_from_p(p_value),
                   "(", fmt_num(std_error, 2), ")")
  ) %>%
  select(treatment, scenario, 조정치) %>%
  pivot_wider(names_from = scenario, values_from = 조정치,
              names_glue = "{scenario} 조정치(SE)")

표3_민감도분석 <- sensitivity_rv_pooled %>%
  left_join(표3_벤치마크, by = "treatment") %>%
  mutate(
    `도서 보유 수` = as.character(`도서 보유 수`),
    `RV(α = .05)`  = fmt_dec(RV, 3),
    `XRV(α = .05)` = fmt_dec(XRV, 3)
  ) %>%
  select(
    `도서 보유 수`,
    `RV(α = .05)`,
    `XRV(α = .05)`,
    ends_with("조정치(SE)")
  )


show_table(
  표3_민감도분석,
  "<표 3> 미관측 교란요인에 대한 도서 보유 수 효과의 민감도 분석",
  "table3.csv"
)

show_note(
  sprintf("주. 기준집단은 '%s'이다. ", book_ref),

  sprintf("RV는 미관측 교란요인이 원인변수와 결과변수의 잔여 변량을 같은 정도로 설명한다고 ",
          "가정할 때, 유의수준 %s에서 추정효과의 통계적 유의성을 제거하는 데 필요한 최소 ",
          "부분설명력이다. ", fmt_dec(SENSE_ALPHA, 2)),

  "극단적 RV는 미관측 교란요인이 결과변수의 잔여 변량을 극단적으로 설명할 수 있다고 가정할 때, ",
  "통계적 유의성을 제거하기 위해 원인변수의 잔여 변량을 설명해야 하는 최소 부분설명력이다. ",
  "RV와 극단적 RV는 Cinelli와 Hazlett(2020)의 부분설명력 모수화에 기초하여 산출하였다. ",
  "미관측 교란편의의 크기는 최종가중치를 적용한 가중 선형회귀의 잔차 기하량을 이용하여 환산하였으며, ",
  "통계적 유의성 판단에는 각 PV의 Fay BRR 분산을 Rubin의 결합규칙으로 통합하여 얻은 표준오차와 ",
  "Barnard-Rubin 자유도를 적용하였다. ",

  "1배, 2배 및 3배 시나리오는 가상의 미관측 교란요인이 비교기준 공변량인 ",
  "'부모 교육수준: 대졸 이상'보다 원인변수 및 결과변수와 각각 1배, 2배 및 ",
  "3배 강하게 관련되는 경우를 나타낸다. ",

  "이 공변량의 관찰된 부분설명력은 결과변수에 대해 R2 = ",
  fmt_dec(r2_y_benchmark, 3),
  ", 원인변수에 대해 범주별로 ",
  fmt_dec(min(r2_d_benchmark), 3),
  "-",
  fmt_dec(max(r2_d_benchmark), 3),
  "였다. 괄호 안은 조정 표준오차이다. ",
  "*p < .05, **p < .01, ***p < .001."
)


# ==============================================================================
# 10. 재현 정보
# ==============================================================================

cat("\n[재현 정보]\n")
cat("  ", R.version.string, "\n", sep = "")
for (pkg in packages) {
  cat(sprintf("   %-12s %s\n", pkg, as.character(packageVersion(pkg))))
}

# ==============================================================================
# 끝
# ==============================================================================
