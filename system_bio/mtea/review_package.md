# 검수 요청 패키지: Benedicto 2025 한계 보완 분석 — 방향 설정 논리와 결과

(검수자용 독립 문서. 이 파일만 읽어도 이해되도록 작성. 상세 수치는 `extension_report.md` 1–10절.)

## 0. 무엇을 검수받고 싶은가
1. 분석 방향(어떤 한계를 어떤 순서로 공략했는가)이 타당한가
2. 각 단계 방법(통계·데이터 해석)에 오류/과잉 해석이 있는가
3. 최종 결론의 주장 강도가 근거에 맞는가
4. 놓친 한계·대안 분석이 있는가

## 1. 대상 논문과 목적
- 논문: Benedicto, Flobak, Ponce-de-León, Valencia 2025, "Constraint based modeling of drug induced metabolic changes in a cancer cell line", npj Syst Biol Appl 11:111. 데이터 GEO GSE285616, 코드 MTEApy.
- 내용: AGS 위암 세포에 PI3K 억제제(PI103), MEK 억제제(PD0325901), TAK1 억제제(5Z-7-oxozeaenol), 그리고 PI3Ki+MEKi, PI3Ki+TAKi 병용(각 약 GI50의 0.5배)을 24 h 처리 → RNA-seq(조건당 4 batch 반복, 총 24 샘플) → DESeq2 → Human-GEM v1.18의 187개 metabolic task에 대해 TIDE(task 반응에 투영한 log2FC 평균), TIDE-essential(task 필수 유전자 log2FC 평균) 점수 → 유전자 라벨 permutation으로 p → 약물별·병용 "metabolic synergy"(weak/strong, 저자 휴리스틱) 주장.
- 내 목적: 중간고사 발표(5분, 논문 분석 설명형) 준비 중 "논문 한계를 보완할 방법"을 찾고 가능한 만큼 실제로 수행. **결과는 PPT에 넣지 않음(사용자 지시), Q&A 방어 자료.**

## 2. 한계 후보 도출과 우선순위 (방향 설정)
논문 정독 후 도출한 한계 후보 → 공개 데이터·짧은 시간 내 수행 가능성과 발표 가치로 순서화.

| # | 한계 후보 | 근거 | 수행 여부 |
|---|---|---|---|
| L1 | 다중검정 보정 없음 (task 187개를 p<0.025로 컷) | 본문·보충표 | 1, 2절 |
| L2 | "synergy"가 통계적 상호작용 검정이 아니라 휴리스틱 (조합 점수가 단일약보다 더 크면 weak, 임계 근처면 strong) | 방법 절 | 3, 5절 |
| L3 | 유전자 라벨 permutation은 유전자 간 상관을 무시 → p 낙관적 | MTEApy 구현 | 4, 6절 |
| L4 | 설계: 병용 반수 용량, 반수 용량 단일약 arm 없음 → 비선형 용량반응이면 "상승효과"가 아티팩트 | 설계 | 3, 7절 |
| L5 | 전사체 점수가 실제 약효(성장)와 연결되는지 미검증 | 논문 | 7, 8절 |
| L6 | 한 층(mRNA)만 봄 (효소 인산화, 대사체 미반영) | 설계 | 9절 |
| L7 | 강조된 후보 유전자(ASNS, ODC1/SRM, PHGDH 축, CoQ, TXNRD1)의 실제 성장 기여 미확인 | 논문 | 10절 |

수행 순서 논리(각 단계 결과가 다음 의문을 낳음):
1. L1: 보정하면 몇 개 남나? → 대부분 탈락 (방법 자체 문제인지, 검정력 문제인지 모름)
2. L2·L4: 병용이 가법(반수 용량 기대치)보다 큰가를 직접 모델링 → 일부 지지, 그러나 L3 때문에 p가 낙관적일 수 있음
3. L3: 상관 인지형 검정(CAMERA-like) + 음성대조(라벨 셔플), 양성대조(spike-in)로 **내 분석 자체를 교정** → 보정된 검정에서 상호작용은 FDR 5%에서 0개
4. HSA(Highest Single Agent; 병용이 두 단일약보다 모두 큰가): 병용 효과는 PI3Ki 단독과 구분 불가 (다만 해상력 한계 명시)
5. L5: 실제 약효 데이터(Flobak 2015 성장 곡선) 대조 → 성장 수준 시너지는 진짜, 그러나 24 h 전사체 점수가 이를 따라가지 않음
6. 독립 데이터(Tsirvouli 2024 시계열 RNA-seq, 동일 계)로 재현성·용량 확인 → 견고한 task 수준 시너지 없음
7. L6: 인산화체로 mRNA 외 층 확인 → mRNA와 무상관
8. L7: DepMap CRISPR·CCLE로 후보 유전자의 실제 필수성 확인 → 후보는 비필수, 핵산 합성이 필수

## 3. 단계별 방법·결과·판정

| 단계 | 질문 | 방법 | 결과(요약) | 판정 |
|---|---|---|---|---|
| 재현 | 논문 점수를 재현하나 | pydeseq2(~batch+condition, control 기준), MTEApy TIDE/TIDE-essential, 10,000 perm | 내 점수 vs 저자: r = 0.95 / 0.94 / 0.78 / 0.82 / 0.40(MEKi), 내 p가 더 보수적 | 대체로 재현 (MEKi는 약함) |
| 1. BH | 보정 후 남는가 | 내 재현에 BH(양측 p=2×꼬리) | TIDE p<.025 → BH q<.05: 22→0, 16→0, 44→2, 29→3, 34→3. TIDE-essential 전 조건 0 | 대부분 탈락 |
| 2. 저자 p에 BH | 저자 보충표 기준 | BH on Supplementary Data 3, 4 | PI3Ki 78→59, PI3Ki+MEKi 64→27, 단일 TAKi/MEKi 0. synergy 48개 중 23(q<.05), strong 8개는 q 0.127–0.158로 전부 불통과, weak 마진 중앙값 0.08 log2 | PI3Ki 효과는 견고, synergy 주장은 취약 |
| 3. 상호작용 | 반수 용량 기대치보다 큰가 | 대비 AB−0.5A−0.5B, 3,000 perm | PI3Ki+MEKi 76 task(음), BH 41; PI3Ki+TAKi 27, BH 0 | (L3 교정 전) 부분 지지 |
| 4. CAMERA-like | 상관을 고려하면 | Wald z, 잔차 상관 추정, VIF=1+(m−1)ρ | 집합 내 중앙 상관 0.10–0.18 → 분산 1.6–2.6배. p<.05: 독립 가정 44–84개 → 추정 ρ 11–15개, BH 0–2 | 낙관적 p 확인 |
| 5. HSA | 병용 > 두 단일약 | 교집합-합집합(max p) | PI3Ki+MEKi 0, PI3Ki+TAKi 0 (BH). 병용 vs PI3Ki 단독은 구분 불가 | "검출 불가" |
| 6. 자기 검증 | 내 검정이 맞는가 | 대비 대수·BH 수치 검증, 라벨 셔플 null, spike-in | 대수/BH 오차 1e-15. 유전자 셔플 FPR 0.39(p<.05) 대 CAMERA-like 0.037. spike-in 검출력(0.3 log2) 80% / BH 42% | CAMERA-like 보정됨, 다소 보수적 |
| 6b. 보정 재검정 | 상호작용 최종 | CAMERA-like 120 task | PI3Ki+MEKi p<.05 40개, BH 0 (최소 q=0.075, serine/glycine 합성); PI3Ki+TAKi 11, 0 | 시사적, 확립 아님 |
| 7. 실제 약효 | 성장 시너지 실재 | Flobak 2015 S1 성장 곡선 digitize (48 h) | 억제율: 병용 78%/92% vs 단일 최대 52–73%. HSA·Bliss 충족. 반수 용량≈0.7–0.8배 효과 | 성장 시너지 실재 |
| 7b. 용량 가중 | 기대치 w를 바꾸면 | AB−wA−wB, w=0.5…1.0 | PI3Ki+MEKi: w 0.5 → 34개(BH 0), w 0.7 → 33(1), w 0.9 → 5(0) | 신호 약해짐 |
| 7c. 증식 서명 | 전사체가 성장 따라가나 | Tirosh S+G2/M 94 유전자 | PI3Ki +0.31, MEKi −0.42, PI3Ki+TAKi +0.25(성장 억제 최대인데 증가) | 전사체≠성장 순위 |
| 8. 독립 데이터 | 재현·용량 | Tsirvouli 94 샘플, 1–24 h, limma식 | 24 h gene log2FC r 0.54–0.79. 용량 진단 기울기 0.47–0.75(GSE285616) vs 0.76–0.95(전량). task 집합 BH 0–3개, 반복 신호는 keratan sulfate/glycan 분해 상향 | 견고한 task 시너지 없음 |
| 9. 인산화체 | mRNA 외 층 | 18,491 site, 560 효소 | 효소 농축 없음(OR 0.84–1.02), 인산화~mRNA ρ −0.13…+0.12, 인산화 변화 효소 37%는 mRNA 변화 없음 | 별개 층 |
| 10. DepMap | 후보 유전자 실제 필수성 | CRISPR Chronos 24Q4, AGS | ASNS −0.06, PHGDH −0.01, PSPH 0.00, ODC1 −0.15, SRM −0.10, TXNRD1 −0.21 (비필수). RRM1/RRM2/DUT/UMPS/DHODH −1.5…−2.3(필수). task 필수 유전자 29.2% vs 12.7% 의존 | 후보 비필수, 핵산 합성 필수 |

## 4. 자기 정정(과잉 주장 철회) 이력
- "반수 용량 병용 = PI3Ki 단독" → "차이 검출 불가(해상력 task당 약 0.2–0.6 log2)". 오히려 8 task에서 병용이 덜 극단적.
- 상호작용 41 task → 보정 후 FDR 5%에서 0 (gene-shuffle p는 낙관적).
- 저자 BH 카운트는 gene-shuffle p 기반이므로 상한선.
- 라벨 셔플 DESeq2 null은 DEG 0이라 퇴화(점수 검정 불가) → set 수준 null을 따로 구성.
- DE 첫 실행에서 reference level 오류(알파벳순 5Z) → control 기준으로 재실행.
- 1,000 permutation은 BH에 지나치게 엄격 → 10,000으로 재실행.

## 5. 최종 결론과 주장 강도

| 주장 | 근거 강도 |
|---|---|
| 성장 수준 시너지는 실재 (Flobak 2015) | 강 (독립 실험, 단 digitize·단일 반복) |
| PI3K 억제가 견고한 주 동인 | 강 (저자 p 기준 BH 후 59 task — 상한선; 내 재현은 BH 후 2 task. 두 데이터에서 방향 일치, 24 h gene log2FC r=0.60) |
| 병용의 task 수준 대사 시너지는 이 설계로 확립 불가 | 중 (검정력·반수 용량 설계 한계; "없다"가 아니라 "입증 안 됨") |
| 전사체 task 점수는 성장(약효)을 대리하지 못함 | 중 (24 h vs 48 h 비교, digitize) |
| 논문 강조 후보 유전자는 AGS 성장 필수 아님 | 중 (배지 차이 가능성) |
| 성장 관련 대사 축은 핵산 합성 | 중~약 (CRISPR + mRNA 방향, 탐색적) |
| mRNA는 대사 조절의 한 층일 뿐 | 중 (인산화~mRNA 무상관) |

## 6. 내가 아는 약점 (검수자가 공격해줬으면 하는 지점)
1. **CAMERA-like는 limma::camera가 아니라 내 단순 구현.** Wald z + 15 residual df로 추정한 집합 내 잔차 상관 → 24 샘플에서 노이즈. 교정은 null/spike-in으로 확인했지만 실제 camera/GSEA와 직접 비교 안 함.
2. **용량 가정.** 병용 = 각 약 0.5 GI50. 선형(w=0.5) 가정은 성장 곡선 digitize(눈대중)로 w 0.7–0.8을 추정. 선형 가정 + 기대치 설정이 결론을 좌우.
3. **digitize 데이터.** 단일 대표 반복, 오차 약 ±0.1, 48 h 성장 vs 24 h RNA-seq 시점 불일치.
4. **다중 분석의 forking path.** 임계값(|logFC|≥0.5, |mRNA|<0.25, 마진 0.2), 컬럼 선택을 사전 고정하지 않음. 많이 돌렸고 일부는 탐색적.
5. **독립 데이터 용량 추론.** Tsirvouli 용량이 "전량"이라는 근거는 프리프린트 서술과 기울기 진단뿐 (메타데이터에서 직접 읽지 않음). 기울기는 x 노이즈로 하향 편향.
6. **DepMap/CCLE.** 배지(Ham's F12 아님), 기준선 대사체(약 처리 전), 필수성은 약 반응의 필수 조건이 아님.
7. **인산화.** 단백질 양으로 보정 안 됨, 수천 site가 유의(전역 이동 가능), 핵심 효소 site 미측정(정보 없음 ≠ 변화 없음).
8. **검정력.** 3–4 반복, HSA 해상력 task당 0.2–0.6 log2 → "없음"이 아니라 "검출 불가"라는 표현 유지 필요.
9. **저자 보충표 불일치.** 시트 카운트 31/19/78/56/64 vs 본문 30/17/77/55/66, TAKi/MEKi 헤더 뒤바뀜, 48 vs 49 synergy.

## 7. 검수 질문 (구체)
1. 유전자 집합 competitive test에서 잔차 상관을 직접 추정해 VIF를 적용한 방식이 camera의 inter.gene.cor 고정값 대비 타당한가? 더 나은 대안(예: rotation test FRY/ROAST)은?
2. 교집합-합집합(max p) HSA 검정은 점수 대비(log2FC 차)에 보수적으로 맞는가? 동등성 검정(TOST)으로 "차이 없음"을 주장할 수 있는가?
3. 반수 용량 설계에서 기대치(w·A + w·B)를 성장 곡선으로 보정하는 방식이 Bliss/Loewe 관점에서 적절한가?
4. 24 h 전사체가 48 h 성장을 따라가지 못한다는 결론은 시점 차이를 감안하면 정당한가?
5. 독립 시계열 데이터에서 기울기 진단으로 용량을 추론한 것을 과잉 해석하지 않았는가?
6. DepMap 필수성을 약 반응 해석에 쓰는 것이 합리적인가, 어디까지 주장 가능한가?
7. 위 결론 표(5절)의 강도 부여가 과하거나 부족한 곳은?
8. 내가 놓친 한계 후보가 있는가?

## 8. 재현 정보
- 분석 코드: `repro.py`, `inter.py`, `camera_like.py`, `taskgenes.py`, `hsa_analysis.py`, `validate_checks.py`, `validate_interaction.py`, `nullde.py`, `dose_weight_scan.py`, `prolif_check.py`, `tsirvouli_analysis.py`, `phospho_analysis.py`, `depmap_stream.py`, `external_check.py`
- 데이터: GEO GSE285616; 저자 보충 Supplementary Data 3/4/6 (`si/`); Flobak 2015 PLoS Comput Biol S1 Text; Tsirvouli 2024 repo (github.com/Eirinits/AGS_synergies); DepMap 24Q4 (figshare 공식 릴리스 미러); CCLE metabolomics (Li et al. 2019)
- 환경: Python (pydeseq2, MTEApy 로컬 설치, scipy), Windows

## 9. 검수 라운드 1 반영 (GPT 검수 → 직접 검증 → 수정)
상세 수치: `extension_report.md` 11절. 각 지적을 파일·원문 소스로 직접 확인함.

### 지적이 맞아 수정한 것
| 항목 | 수정 |
|---|---|
| CAMERA-like 자유도 | limma `camera` 소스 확인: 추정 ρ → df = min(잔차 df, G−2). df=15로 재실행: 용량 보정 PI3Ki+MEKi 상호작용 p<.05 40→33개, BH 0 유지. 결론 불변, 숫자만 약간 낙관적이었음 |
| "TIDE 41→0" | 비교 단위가 다름. 187 task 중 56개 <5 유전자, 11개 >200 → 120개 검정(고유 집합 116). "별도 120 task gene-union 검정에서 BH 5% 통과 없음"으로 표현 변경 |
| 해상력 0.2–0.6 | 무작위 유전자 집합 분위수일 뿐 복제 기반 검출한계 아님 → "작은 효과 검정력 미평가" |
| TIDE(AB−A) ≠ TIDE(AB)−TIDE(A) | GPR min/absmax·마스킹 비선형. 유전자 수준 대수 검증만 했음 |
| w-스캔 | Bliss/Loewe 아님, 민감도 분석. 일부 w에서 BH 1개 남아 "사라짐" 과장 |
| 전사체≠성장 | 증식 서명 하나와 digitize 값 비교였음 → "연결 미검증" |
| Tsirvouli 용량 | Methods에 PI .70 / 5Z .50 / PD .035 µM 명시. 기울기로 추론 불필요 |
| Tsirvouli 척도 | 공식 스크립트 확인: VSN(glog2) 출력. 내 "log2FC"는 VSN 척도 차이(저발현 유전자는 0쪽으로 압축; 유전자 23.5%가 평균<7). 공식 limma logFC 미사용 |
| CoQ | COQ2 −0.530, COQ3 −0.502, COQ5 −0.514, COQ4 −0.949, PDSS1 −0.548 → "경계~의존"이지 비필수 아님 |
| 23/48, strong q | 병용 vs control의 BH이지 synergy 검정 아님 |
| 37% | 시점 불일치. 동시점(8 h) mRNA로 재계산 44% |
| 저자 BH "상한선" | 보장 없음 → "낙관적일 가능성" |

### 새로 확인된 것
- **null 범위에 따라 1절 결과가 달라짐** (10,000 perm, BH 5% 통과): PI3Ki 2/2, MEKi 0/3, TAKi 0/0, PI3Ki+TAKi 3/5, PI3Ki+MEKi 3/20 (전체 유전자 null / 모델 유전자 null[MTEApy 기본]). 저자: 59/0/0/3/27.
- 그래서 PI3Ki+MEKi 격차(20 vs 27)는 작고, **PI3Ki 격차(2 vs 59)만 큼**. 점수는 r=0.95로 일치하나 내 null 폭이 저자보다 넓음(암시 SD 0.26 vs 0.19). shrinkage 제거·마스킹 제거·null 범위 변경 모두 격차를 못 닫음 → **원인 미해결**.

### 수정된 결론 강도
| 주장 | 강도 | 문구 |
|---|---|---|
| 성장 수준 synergy | 선행 실험 근거 중–강 | "Flobak 2015(CI 분석)에서 보고됨. 내 digitize는 보조이며 정밀도 미검증" |
| PI3K 주 동인 | 중 (범위 축소) | "PI103 포함 조건에서 공통 전사 반응 뚜렷. PI103은 mTOR/DNA-PK도 억제하므로 순수 PI3K 인과성은 미확인" |
| task 수준 synergy 확립 불가 | 중 (유지) | "현재 설계·분석으로 확립 못함. 부재를 증명한 것은 아님" |
| 전사체 task 점수가 성장 대리 | 약 | "연결 미검증" |
| 후보 유전자 | 중 | "강한 기준선 의존성 없음. CoQ는 경계. 약물 조건 의존성은 별개" |
| 핵산 합성 축 | 약 (탐색) | "강한 기준선 의존성 + 일부 mRNA 감소. synergy 매개 여부 미확인" |
| mRNA 한 층 | 일반 명제는 타당, 본 분석은 보조 | "인산화~mRNA 순위 상관 낮음. 층간 독립성·flux 변화 입증 아님" |

### 허락 후 추가 수행 (extension_report 12–13절)
1. **공식 limma logFC로 Tsirvouli 재확인**: 공식 vs VSN 차이 r=0.82–0.90, VSN은 저발현 유전자(평균<7)에서 기울기 0.26–0.37로 압축, 발현 좋은 유전자는 0.86–0.90. ASNS, 증식 서명, 핵산 합성 군 결론 방향 불변(크기는 VSN이 0–15% 작게 나옴). **단, 용량 기울기 진단은 반량/전량을 구분 못함**: 전량 데이터 24 h 기울기 0.46–0.86 vs GSE285616 0.47–0.75. "반량 증거"라는 이전 서술 철회; 반량은 논문 Methods와 Flobak 2015에 근거.
2. **limma-voom camera / fry (R 4.6.1)**: 같은 count, 같은 120 task 집합, 같은 대비.
   - 추정 ρ camera는 내 CAMERA-like와 거의 일치 (용량 보정 PI3Ki+MEKi 30 vs 33개, vs PI3Ki 8 vs 8 전부 상향, vs MEKi 79 vs 85) → 내 교정이 맞았음.
   - **limma 기본값 ρ=0.01은 이 데이터에서 반보수적**: 라벨 셔플 FPR 0.28–0.35, 셔플당 BH 24–29개 거짓 양성. 기본 camera로는 상호작용 81/120 task가 "유의"하게 나옴 → 신뢰 불가.
   - fry(FPR 0.05–0.07, 거의 보정됨): PI3Ki+MEKi는 w=0.7–0.8에서도 33–36 task가 반수 용량 기대치에서 벗어남. 하지만 경쟁검정(추정 ρ)은 모든 w에서 BH 0 → **편차는 광범위한 전사 비가법성이며 task 특이적이지 않음.** PI3Ki+TAKi는 거의 없음(1–4).
   - 저자 PI3Ki+MEKi synergy 40개 중 39개가 기본 camera BH 유의이지만 전체 집합도 81/120이 유의 → 농축 거의 없음.
3. 새 결론 문구: "PI3Ki+MEKi 반수 용량 병용은 선형 가법 기대치와 다른 광범위한 전사 비가법성을 보인다(fry). 그러나 특정 task의 시너지로 국한되지 않는다(경쟁검정 0). 이 비가법성이 비선형 용량반응에서 오는지 진짜 상호작용인지는 반량 단일약 arm 없이 구분 못한다."

### 아직 못 하는 것
- 반량 단일약 RNA, 동시점 성장: 공개 자료 없음

## 10. 도구(MTEApy/TIDE) 신뢰도 검증과 췌장암 적용 (extension_report 14절)
| 검사 | 결과 | 판정 |
|---|---|---|
| 양성 대조: 스타틴 → 콜레스테롤/메발론산/FPP 합성 과제 상향 예상 (사전 지정 X138, X140, X141) | 20 µM PANC-1 세 스타틴 모두 1–4위(p≤0.009). MiaPaCa-2는 메발론산 2–3위, 나머지 11–47위. 2 µM은 약함 | 강한 전사 반응에서만 작동 |
| 방향 | 효소가 억제되는데 점수는 상향 | mRNA 보상 반응이지 flux 방향 아님 |
| 독립 실험 재현 | 같은 약 ρ 0.60–0.72, 다른 약 평균 0.53. TIDE-essential 더 약함 | 재현되나 약 특이성 약함 |
| GPR 투영의 부가가치 | 단순 평균과 r 0.88–0.93 | 거의 없음 |
| p-값 | 유전자 라벨 셔플 FPR 약 0.39 | 사용 불가, 점수·순위만 |
| 세포주 간 재현 | PDAC 4주 평균 ρ 0.20 (PANC-1 거의 0) | 세포주 의존 큼 |
| 교차암 약 특이성 | PDAC trametinib vs AGS MEKi 0.42–0.51, vs PI3Ki/TAKi 0.10–0.34 | MEKi 패턴 공유 |
PDAC trametinib 24 h 합의 과제(3/4주 이상): 핵산(dCTP/dTTP/dATP/CTP/UTP/GTP) 합성·복합체 I·spermidine 하향, 콜레스테롤·FPP 합성 상향. AGS MEKi 두 데이터와 21/24 방향 일치. DepMap 췌장 라인에서 핵산 과제 100%, 복합체 I 96%, FPP 89% 라인이 의존. 단 증식 신호(S기 서명에 TYMS/RRM1/RRM2 포함)와 분리 불가, 병용·flux·대사체 데이터 없음 → 가설 수준.

## 11. flux 방향 검증: E-Flux가 실제 성장 억제를 재현하나 (extension_report 15절)
- 방법: Human-GEM + Ham's F12 유사 배지, 발현으로 반응 용량을 제한한 FBA(E-Flux), GSE285616 24 h 샘플별로 biomass 예측 → 대조군 대비 비율 → Flobak 2015 48 h 성장 억제(47–92%)와 비교.
- 결과: 18개 설정 모두 예측 억제율이 −19%~+5% (대부분 성장 증가 또는 무변화). 실제는 47–92%.
- 구현 점검: RRM1/RRM2, TYMS, DHODH, PHGDH 축, ASNS, ODC1, HMGCR, GAPDH를 90% 낮춰도 성장 변화 0.0% → 이 설정에서 모델이 해당 효소에 둔감함.
- 결과 의존성: 용량 스케일과 미검출 유전자 floor라는 자유 파라미터 2개에 의존, 데이터로 식별 불가. lipid/cofactor pool은 배지로 생산 불가여서 제외.
- 판정: mRNA에서 flux 방향·성장을 추정하는 것은 이 파이프라인으로는 근거 없음. 단 단일 generic 모델·E-Flux·미보정 배지라는 한계가 있어 "모든 flux 방법이 안 된다"는 뜻은 아님.

## 12. 성장 수준 시너지 판정: 표준 모델로 재계산 (extension_report 16절)
- 데이터: Flobak 2015 S1 Text의 21개 병용 성장 곡선(48 h)과 단일약 용량-반응 곡선을 픽셀 단위로 디지타이즈. 용량 행렬이 없어 MuSyC/BRAID는 불가.
- 지표: HSA 마진(저자 기준), Bliss 초과, Loewe 조합지수(CI).

| 조합 | f(A)/f(B) | f(AB) | HSA 마진 | Bliss 기대 | Bliss 초과 | CI |
|---|---|---|---|---|---|---|
| PI3Ki+MEKi | 0.52/0.48 | 0.77 | +0.25 | 0.59 | +0.18 | 0.26 |
| PI3Ki+TAKi | 0.76/0.61 | 0.96 | +0.20 | 0.71 | +0.25 | 0.25 |

- 21개 조합 중 저자가 시너지로 부른 4개와 아닌 17개를 HSA·Bliss·CI 모두 분리 (AUC 1.00, CI는 유효한 10개에서 완전 분리).
- 같은 약의 실험 간 변동 sd 0.08–0.11. 비시너지 17개의 Bliss 초과 −0.13±0.18 대비 PI3Ki+MEKi는 약 1.8 sd, PI3Ki+TAKi는 약 2.2 sd.
- 민감도: 반량/전량 효과비 r이 1에 가까워지면(용량 반응이 매우 평평) 초과가 사라짐. 측정 곡선은 r≈0.66–0.8.
- 판정: 성장 수준 시너지는 두 조합 모두 지지됨(중간 크기, PI3Ki+TAKi가 더 확실). 이 시너지는 전사체 task 점수로는 확인되지 않음.

