# Benedicto 2025 보완 분석: 방향은 타당하고 결론은 좁혀야 한다

**분석 방향은 타당하나, 현재 결과로 task 시너지를 확정할 수 없다.** 병용 대 control의 BH, TIDE 대비 점수, gene-union 경쟁검정, 성장 synergy를 같은 증거로 연결한 점이 핵심 오류다.

대상은 AGS 한 세포주, 24 h RNA-seq, 6조건×4반복이다. 단일약은 GI50, 병용은 각 약의 0.5 GI50다. TIDE는 유전자 발현차를 GPR로 반응에 투영한 평균 점수다. 대사 flux나 성장률을 직접 측정하지 않는다. ([Benedicto Methods](https://www.nature.com/articles/s41540-025-00586-y), [TIDE 구현](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/tide.py:196))

검수는 원 논문·공식 코드, 정적 코드, 저장 출력, 일부 표·행 확인에 한정했다. 전체 재실행·코드 수정은 하지 않았다. 수치 확인과 검수 판단을 구분한 Q&A 자료이며 PPT에는 반영하지 않는다.

## 1. 잔차 상관 추정과 VIF 적용은 CAMERA로 인정할 수 있는가?

**원리는 타당하다. 현재 구현을 공식 CAMERA 또는 교정된 TIDE 검정으로 부를 수는 없다.** 평균 잔차 상관과 `VIF=1+(m−1)ρ`는 CAMERA 원리에 맞는다. `inter.gene.cor=0.01`은 보편적 진실값이 아니다. 공식 도구도 `inter.gene.cor=NA`로 집합별 상관을 추정한다. 직접 추정했다는 사실 자체가 오류는 아니다. ([CAMERA 원 논문](https://pmc.ncbi.nlm.nih.gov/articles/PMC3458527/), [limma 매뉴얼](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf))

**자유도 처리에는 차이가 있다.** 공식 CAMERA는 고정 ρ에 `G−2`, 추정 ρ에 `min(residual df,G−2)`를 쓴다. 현재 설계의 잔차 자유도는 15다. 로컬 코드는 ρ를 추정하면서 유전자 수 기반 자유도를 쓴다. `G−2`가 항상 틀린 것은 아니다. 추정 ρ의 불확실성을 고정 ρ처럼 처리한 점이 문제다. 또한 DESeq2 Wald z, 배경 유전자 분산, 별도 무가중 log-count OLS 잔차를 결합한다. 공식 CAMERA의 moderated 통계량과 pooled variance 절차와 다르다. ([공식 CAMERA 소스](https://raw.githubusercontent.com/bioconductor-source/limma/master/R/geneset-camera.R), [로컬 계산](D:/Cander/School/3-2/system_bio/mtea/camera_like.py:34), [최종 계산](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:35))

검정 질문도 바뀌었다. CAMERA-like는 집합의 평균 z가 나머지 유전자보다 다른지 검정한다. TIDE의 반응별 GPR 점수와 다르다. 최종 120 task는 유전자 수 5–200 필터를 거친 gene union이다. GPR의 AND/OR 논리와 반응별 가중을 보존하지 않는다. 따라서 “TIDE 상호작용 41개를 교정하니 0개” 대신 **“별도 120 task gene-union 경쟁검정에서 BH 5% 통과 집합이 없었다”**고 써야 한다. 최소 q=.075는 저장 보고서의 기록이며 이번에 독립 재계산하지 않았다. ([집합 생성](D:/Cander/School/3-2/system_bio/mtea/taskgenes.py:7), [크기 필터](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:18), [출력 기록](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:98))

null·spike-in은 유용한 내부 점검이다. 다만 FPR .387/.037과 spike 0.3의 raw/BH 검출 .80/.42는 unmasked OLS t→z 기반이다. 최종 DESeq2 Wald z나 masked TIDE/HSA 파이프라인의 오류율·검정력을 검증한 값은 아니다. 실제 DESeq2 null은 한 셔플에서 DEG 0으로 퇴화했다. **“OLS 합성 null에서 해당 오류율을 관찰했다”**가 정확하다. ([검증 계산](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:33), [로그](D:/Cander/School/3-2/system_bio/mtea/validate_checks.log:16), [full-pipeline null](D:/Cander/School/3-2/system_bio/mtea/nullde.py:6))

권고: 동일 universe·집합·설계·대비로 limma-voom CAMERA를 비교한다. FRY/ROAST는 집합 자체가 변했는지 묻는 self-contained 검정이다. 상관을 보존하는 보완 분석으로 적절하지만 CAMERA와 귀무가설이 같지는 않다. ([limma 매뉴얼](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf))

## 2. max p HSA는 맞는가? TOST로 차이 없음을 주장할 수 있는가?

**사전 고정한 방향의 단측 p라면 `max(pA,pB)`가 맞다.** 방향을 s로 두고 `δA=s(AB−A)`, `δB=s(AB−B)`라 정의한다. 대립가설은 두 차이가 모두 양수다. 구성 검정 하나라도 실패하면 HSA도 실패한다. 이 IUT는 두 p의 독립성을 요구하지 않는다. 공유 병용·control의 공분산은 각 대비의 SE를 모델에서 계산할 때 반영해야 한다. 현재 유전자 대비는 DESeq2 계수 대비로 계산하며, 단순 독립 SE 합산은 아니다. ([IUT 원 논문](https://doi.org/10.1214/ss/1032280304), [계수 대비](D:/Cander/School/3-2/system_bio/mtea/inter.py:35))

현재 구현은 max 양측 p와 사후 sign 필터를 쓴다. 대칭 null과 올바른 방향 조건에서는 보수적일 수 있으나 정식 directional HSA로 인증할 근거는 부족하다. 반대 방향인 mixed task도 BH 후 제거하므로 HSA 가설에 맞는 p를 먼저 구성해야 한다. 더 큰 문제는 **`TIDE(AB−A)≠TIDE(AB)−TIDE(A)`**라는 점이다. min/absmax GPR과 대비별 유의성 마스킹은 비선형이다. 유전자 대비 대수 오차 1e−15는 task 점수 차나 그 p의 타당성을 보장하지 않는다. ([HSA 코드](D:/Cander/School/3-2/system_bio/mtea/hsa_analysis.py:17), [마스킹·GPR](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/utils.py:112))

TOST는 적용 가능하지만 현재 자료로 동등성을 주장할 수 없다. 생물학적으로 정당한 한계 ±Δ를 결과 확인 전에 정해야 한다. α=.05에서 해당 효과 차의 90% CI가 구간 안에 들어야 한다. 현재 비유의 경쟁검정과 random-gene null만으로 이 CI를 얻을 수 없다. **0.2–0.6 “해상력”은 masked LFC random-gene 평균의 분위수다. 생물학적 반복 기반 최소 검출효과나 TOST 한계가 아니다.** 문구를 삭제하고 “추가 효과 미검출; 작은 효과의 검정력 미평가”로 바꾼다. ([TOST 방법 설명](https://doi.org/10.1177/1948550617697177), [해상력 계산](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:126))

## 3. 성장 곡선으로 w를 보정하면 Bliss/Loewe 기대치인가?

**아니다. 전사체 가법 대비의 가중치 민감도 분석이다.** Bliss는 같은 endpoint·투여량에서 `E_A+E_B−E_AE_B`를 기대 억제율로 쓴다. Loewe는 단일약 용량반응에 따른 용량 등가성을 쓴다. 성장 억제율의 반량/전량 비율을 mRNA log2FC에 곱하면 둘 중 어느 모델도 되지 않는다. 성장의 w와 전사체의 w가 같다는 연결모형 또는 반량 단일약 RNA가 필요하다. ([SynergyFinder 2.0 모델 정의](https://academic.oup.com/nar/article/48/W1/W488/5815821), [가중 대비 코드](D:/Cander/School/3-2/system_bio/mtea/dose_weight_scan.py:30))

Flobak 2015는 48 h 임피던스 cell index와 Chou–Talalay CI로 성장 synergy를 보고했다. 선행 근거는 사용자 디지타이즈 한 곡선만이 아니다. 다만 78%/92%와 w≈.7–.8의 정밀도는 원좌표·반복별 오차 전파가 없어 검증되지 않았다. 성장 synergy 근거와 특정 metabolic task의 기전 근거를 분리해야 한다. ([Flobak 원 논문](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004426), [대표 반복 설명](D:/Cander/School/3-2/system_bio/mtea/si/flobak_S1.txt:782))

반량 병용을 전량 단일약과 비교하는 HSA는 투여량이 일치하지 않는다. 같은 endpoint에서 단조 용량반응이 성립하면 전량 단일약을 넘는 결과는 보수적 근거다. 넘지 못했다면 반량 단일약 HSA를 배제할 수 없다. task 변화의 단조성도 미확인이다. w 스캔은 비단조이고 일부 w에서 BH 1개가 남으므로 “보정하면 사라짐”도 과하다. **“선택한 w에 결과가 민감했다”**가 맞다. ([용량·결과표](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:115))

## 4. 24 h 전사체가 48 h 성장을 대리하지 못한다는 결론은 정당한가?

**“현재 분석에서 두 지표의 연결을 검증하지 못했다”까지다.** 실제 비교 대상은 Tirosh S/G2M 유전자 평균과 선행 성장값이다. metabolic task 점수의 예측 성능을 평가하지 않았다. RNA와 성장은 동일 샘플·시점·반복으로 짝지어지지 않았고, 서로 다른 성장실험의 PI103 억제율 52%와 73%도 한 순위에 들어갔다. ([증식 서명 계산](D:/Cander/School/3-2/system_bio/mtea/prolif_check.py:4), [고정 입력 성장값](D:/Cander/School/3-2/system_bio/mtea/prolif_check.py:37), [비교 기록](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:124))

24 h 발현은 세포주기 정지 위치·세포 구성·적응을 반영할 수 있다. 48 h 누적 cell index와 순위가 달라도 대리성 전반의 실패를 뜻하지 않는다. 임피던스는 직접 세포수·세포사멸 측정과도 구별해야 한다. 권장 문장: **“24 h 증식 서명의 조건 순위는 선행 48 h 성장실험과 일치하지 않았다. task 점수의 약효 대리성은 미검증이다.”** ([Flobak Methods](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004426))

## 5. 독립 시계열의 기울기로 용량을 추론한 것은 과한가?

**용량 추론은 불필요하며, 기울기의 약리학적 해석은 과하다.** Tsirvouli Methods는 PI .70 µM, 5Z .50 µM, PD .035 µM를 단일·병용에 투여했다고 명시한다. “전량”은 이 기록으로 판단한다. 기울기는 측정오차, 선택 편향, 공유 DMSO, 시간반응, 상호작용의 영향도 받는다. partner 단독 반응이 작다는 이유로 병용에서 partner 효과가 없다고 가정할 수 없다. ([Tsirvouli Methods](https://www.biorxiv.org/content/10.1101/2024.03.12.584561v1.full), [기울기 선택·회귀](D:/Cander/School/3-2/system_bio/mtea/tsirvouli_analysis.py:60))

**입력 척도 오류를 먼저 수정해야 한다.** 내려받은 정규화 행렬은 VSN 변환값이다. VSN은 큰 값에서 log2에 접근하지만 일반 log2와 같지 않다. 현재 평균차·기울기·log2 기반 임계값을 정확한 log2FC 분석으로 해석하면 안 된다. 공식 별도 DE는 TMM–voom–limma다. 공식 logFC 출력으로 핵심 수치를 확인하거나 “VSN 변환 발현차”로 라벨을 바꾼다. ([공식 생성·DE 코드](https://github.com/Eirinits/AGS_synergies/blob/main/Transcriptomics/scr/tr_limma_per_Timepoint.R), [VSN 공식 설명](https://bioconductor.org/packages/release/bioc/vignettes/vsn/inst/doc/A-vsn.html), [공식 DE 출력](https://github.com/Eirinits/AGS_synergies/blob/main/Transcriptomics/Output/20230904_DEGs_perContrast_all.csv))

같은 연구팀의 별도 AGS 실험이며 병용 용량도 다르다. **“별도 전량 투여 시계열에서 방향 일치를 탐색했다”**가 적절하다. 외부 기관의 dose-matched 검증은 아니다. 94샘플은 처리 89개와 untreated 5개이며 PD 1 h는 n=2다. 반복 라벨이 batch인지 확인 후 모델 반영 여부를 정해야 한다. 생물시료 ID 수준의 완전한 비재사용 확인은 남아 있다. ([공식 metadata](https://github.com/Eirinits/AGS_synergies/blob/main/Transcriptomics/Data/sample_info_to_use.tsv), [로컬 모델](D:/Cander/School/3-2/system_bio/mtea/tsirvouli_analysis.py:24))

## 6. DepMap 필수성은 약물 반응 해석에 어디까지 유효한가?

**기준선 성장 의존성의 참고 근거다. 약물 조건의 기전 반증은 아니다.** Chronos는 knockout에 따른 상대 성장효과다. 비필수 중앙값 0, 공통 필수 중앙값 −1로 정규화하며 −.5는 본 분석의 실용 cutoff다. ASNS·PHGDH·PSPH·ODC1·SRM·TXNRD1의 저장 점수 −.21…0은 재확인됐다. 다만 완전 knockout, 부분 mRNA 감소, 약물 처리 중 적응 의존성은 같은 개입이 아니다. ([Chronos 공식 설명](https://raw.githubusercontent.com/broadinstitute/chronos/master/README.md), [AGS 저장 행](D:/Cander/School/3-2/system_bio/mtea/external/crispr_gastric_rows.tsv), [cutoff](D:/Cander/School/3-2/system_bio/mtea/external_check.py:14))

**CoQ 전체를 비필수로 묶으면 오류다.** COQ2 −.530, COQ3 −.502, COQ5 −.514는 자체 cutoff를 통과한다. 경계 수준 의존성으로 따로 쓴다. RRM1/RRM2/DUT/UMPS/DHODH −1.5…−2.3은 강한 기준선 의존성이지만 병용 성장 억제의 매개기전이라는 인과 증거는 아니다. ([CoQ 수치·결론](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:162))

Model_24Q4의 AGS 배지는 F-12K+10% FBS이며 원 연구는 Ham’s F12+5% FCS다. 실제 해당 CRISPR screen 배지는 메타데이터만으로 확정하지 못했다. CCLE 대사체도 약물 처리 전 자료다. 권장 문장: **“지정한 후보는 해당 DepMap 조건에서 강한 기준선 성장 의존성을 보이지 않았다. CoQ 일부는 경계 수준 의존성을 보였다.”** ([AGS 메타데이터](D:/Cander/School/3-2/system_bio/mtea/external/Model_24Q4.csv), [CCLE 원 논문](https://www.nature.com/articles/s41591-019-0404-8))

## 7. 최종 결론 7개의 강도는 어떻게 고쳐야 하는가?

먼저 **23/48은 synergy 검정 통과 수가 아니다.** 저자 synergy-call 중 23개에서 병용 대 control p의 BH가 통과했다는 뜻이다. strong 8개의 q=.127–.158도 같은 병용 대 control 검정이다. “단일약 비유의·병용 유의”는 두 조건 차이가 유의하다는 근거가 아니다. 저자 p의 BH 계산은 맞아도 잘못된 synergy 정의를 해결하지 않는다. ([계산 대상](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:27), [Benedicto synergy 정의](https://www.nature.com/articles/s41540-025-00586-y))

**PI103의 BH 59 대 재현 2, MEKi 점수 r=.40은 미해결 재현 격차다.** “점수 방향·순위는 대체로 유사하나 유의성은 미재현”으로 적는다. 두 분석 모두 10,000 permutation을 사용했으므로 횟수 차이가 설명은 아니다. 저자 BH 카운트는 엄밀한 상한선도 아니다. 검정 통계량·배경이 달라 보정 후 개수가 반드시 줄어든다는 보장은 없다. ([재현 코드](D:/Cander/School/3-2/system_bio/mtea/repro.py:21), [재현·BH 기록](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:19))

| 기존 결론 | 권장 강도 | Q&A용 수정 문장 |
|---|---|---|
| 성장 수준 synergy 실재 | 선행 실험 중–강; 새 디지타이즈는 보조 | “선행 성장실험에서 synergy가 보고됐다. 본 디지타이즈의 정밀도·신뢰구간은 미검증이다.” |
| PI3K 억제가 견고한 주 동인 | 중 | “PI103 포함 조건에서 공통 전사 반응이 뚜렷하다. task 유의성 재현과 PI3K 선택적 인과성은 미확인이다.” |
| task 대사 synergy 확립 불가 | 중; 결론 타당 | “현재 설계와 분석으로 task 수준 synergy를 확립하지 못했다. 효과 부재를 증명하지 않았다.” |
| task 점수는 성장을 대리 못함 | 약; 범위 축소 | “24 h 증식 서명과 선행 48 h 성장 순위가 불일치했다. task 점수의 대리성은 미검증이다.” |
| 강조 후보는 AGS 성장 비필수 | 지정 후보의 기준선 결과 중 | “해당 DepMap 조건에서 강한 기준선 의존성이 없었다. CoQ 예외와 약물 조건 의존성이 남는다.” |
| 성장 관련 대사 축은 핵산 합성 | 약; 탐색 후보 | “핵산 합성 유전자는 강한 기준선 의존성과 일부 전사 감소를 보였다. synergy 매개 여부는 미확인이다.” |
| mRNA는 대사 조절의 한 층 | 일반 명제 타당; 본 분석은 보조 | “관측 phosphosite와 mRNA 변화의 단조 상관이 낮았다. flux 변화나 두 층의 독립성을 증명하지 않았다.” |

PI103은 mTOR/TORC1·DNA-PK 활성도 억제한다. 따라서 관측 약물 반응을 순수 PI3K 경로의 인과성으로 옮기면 안 된다. ([PI103 약리학 원 논문](https://aacrjournals.org/cancerres/article/67/12/5840/533102/Pharmacologic-Characterization-of-a-Potent))

## 8. 놓친 한계는 무엇이며, 다음 분석은 무엇부터 해야 하는가?

**새 분석보다 재현과 검정 대상의 일치가 먼저다.** 저자·재현 간 차이를 gene ID/중복 처리, 필터, shrinkage, masking, Human-GEM 매핑, null universe, p 계산 단계로 나눠야 한다. 모델 유전자와 전체 발현 유전자 배경을 바꾸면 null도 바뀐다. 원 permutation 구현은 p=0을 허용한다. 표본추출 방식에 맞는 +1 보정과 Monte Carlo 불확실성을 점검해야 한다. 1,000에서 10,000으로 늘리는 것은 정밀도 개선이며, “1,000은 항상 BH에 과도하게 엄격”은 성립하지 않는다. ([null 배경 선택](D:/Cander/School/3-2/system_bio/mtea/repro.py:49), [p 계산](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/utils.py:238), [permutation p 원 논문](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf))

task 수도 독립 경로 수가 아니다. 187→120에서 67개가 제외됐고 120개 task 중 고유 gene set은 116개다. 중복 자체로 BH가 무효인 것은 아니나 여러 task를 독립 재현처럼 세면 안 된다. 조건·시간·방법·w를 탐색한 뒤 가장 좋은 결과를 고르면 개별 BH가 전체 탐색을 통제하지 않는다. task별 최소 Chronos 순위도 큰 집합에 유리하다. ([task 집합](D:/Cander/School/3-2/system_bio/mtea/out/task_genes.json), [task 의존성 계산](D:/Cander/School/3-2/system_bio/mtea/external_check.py:37))

인산화의 37%는 8 h site 변화와 24 h mRNA |차이|<.25를 연결한 값이다. 동시점 “mRNA 변화 없음” 검정이 아니다. 동일 유전자의 site 평균은 활성·억제 효과를 상쇄할 수 있다. site 단위 농축에는 단백질별 중복관측이 들어가며 total protein 보정도 없다. OR≈1과 낮은 ρ는 효소 조절 부재·층간 독립성의 증거가 아니다. ([시간·site 계산](D:/Cander/School/3-2/system_bio/mtea/phospho_analysis.py:22), [37% 정의](D:/Cander/School/3-2/system_bio/mtea/phospho_analysis.py:52))

마지막으로 모델의 task-essential은 실제 성장 필수성과 다르다. pFBA 반응 선택, 대안 flux 해, 배지·영양물 공급에 따라 모델 해석이 달라진다. mRNA 감소는 flux 감소나 성장 억제 원인을 자동으로 뜻하지 않는다. AGS 한 세포주 일반화도 제한된다. ([Benedicto TIDE-essential·대안 해](https://www.nature.com/articles/s41540-025-00586-y))

가장 경제적인 후속 순서는 세 가지다. 첫째, estimand·척도·universe·매핑을 고정하고 59 대 2의 격차를 분해한다. 둘째, batch를 반영한 표준 voom–CAMERA와 FRY를 같은 대비로 비교하고 실제 사용 파이프라인의 null을 검증한다. 셋째, Tsirvouli의 VSN/logFC 라벨과 용량 기록을 바로잡고 탐색 결과의 범위를 제한한다. 반량 단일약 RNA, 동시점 성장, 약물 조건 CRISPR/rescue, 실제 대사 flux는 새 실험이다. 공개 자료 재분석으로 누락 arm을 복원할 수 없다.

**발표 Q&A용 최종 문장:** “다중검정과 유전자 상관을 점검하는 보완 방향은 타당했지만, 원 연구의 task 수준 synergy는 현재 설계와 분석으로 확립되지 않았습니다. 추가 경쟁검정의 BH 유의 신호 부재는 synergy 부재나 PI103 단독과의 동등성을 뜻하지 않습니다. 선행 성장 synergy와 전사체 기반 task 변화의 연결은 반량 단일약과 동시점 기능 측정으로 검증해야 합니다.”
