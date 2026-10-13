# 반수 용량·성장 약효·Tsirvouli 시계열 검수

## Q3. 성장 곡선으로 w·A+w·B를 보정하면 Bliss/Loewe 기대치인가?

### Takeaway
아니다. 성장 억제율의 절반/전량 비율을 24 h mRNA log2FC에 이식한 가중치는 약리학적 기대치로 검증되지 않았다. 현재 분석은 임의 전사체 가법 대비의 가중치 민감도 분석으로만 제시해야 한다.

### Cited Findings
- `dose_weight_scan.py`는 성장곡선 비율 약 0.7–0.8을 선언한 뒤, 각 유전자 발현계수에 AB−wA·A−wB·B 대비를 적용한다. 실제 반수 용량 전사체 arm은 없고, 유전자별·task별 용량반응은 추정하지 않는다. [코드 1–2, 30–40줄](D:/Cander/School/3-2/system_bio/mtea/dose_weight_scan.py:30)
- Bliss는 같은 endpoint와 투여량에서 독립적인 단일약 반응의 곱을, Loewe는 같은 물질을 합쳤을 때의 용량 등가성을, HSA는 같은 투여량에서 가장 큰 단일약 반응을 기준으로 쓴다. [SynergyFinder 2.0 원논문](https://academic.oup.com/nar/article/48/W1/W488/5815821)
- Benedicto는 단일약 GI50 PI103 0.700 µM, PD0325901 0.035 µM, 5Z 0.500 µM, 병용 각 0.5 GI50, 24 h RNA 측정, 조건별 n=4를 사용한다. [Benedicto Methods·Table 2](https://www.nature.com/articles/s41540-025-00586-y); [로컬 PDF 추출 695–721줄](D:/Cander/School/3-2/system_bio/mtea/paper_pdf.txt:695)
- Flobak 2015는 xCELLigence 임피던스 cell index 성장분석과 48 h Chou–Talalay CI를 사용했다. 단일약·병용을 한 생물학적 반복 안에서 각각 3단계 희석으로 비교했고, 병용 최고농도는 각 약 0.5 GI50였다. [Flobak 원논문 Methods](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004426); [S1 Text 782–816줄](D:/Cander/School/3-2/system_bio/mtea/si/flobak_S1.txt:782)
- 원논문은 PI3K–MEK·TAK1–PI3K 등의 성장 synergy를 CI가 0.5보다 훨씬 작았다고 보고한다. 따라서 성장 synergy에 대한 선행 근거는 사용자 digitize의 단일 대표 반복만이 아니다. [Flobak Results, Experimental validation](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004426)
- 보고서의 w 스캔은 단조롭지 않다. MEK 병용은 w=.5/.6/.7에서 nominal 유의 task 34/38/33, BH 0/0/1이며 TAK 병용도 .6–.8에서 BH 1이다. [extension_report 115–120줄](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:115)

### Inferences
- 억제율 E에 대한 Bliss 기대치는 E_A(d_A)+E_B(d_B)−E_A(d_A)E_B(d_B). log2 전사체 A+B는 발현비의 곱을 가정하지만, 성장 억제율과 같은 척도가 아니다. w_growth를 w_RNA로 정당화하려면 별도의 연결모형 또는 반수 용량 RNA 측정이 필요하다. [SynergyFinder 모델 정의](https://academic.oup.com/nar/article/48/W1/W488/5815821); [실제 구현](D:/Cander/School/3-2/system_bio/mtea/dose_weight_scan.py:32)
- 용량이 다른 HSA 비교에서 반수+반수 병용이 전량 단일약보다 강하면, 같은 endpoint·실험에서 단일약 효과의 단조성이 성립할 경우 반수 단일약 HSA보다도 강하다는 보수적 근거가 된다. 반대로 전량 단일약을 못 넘었다고 matched half-dose HSA를 배제할 수는 없다. 전사체 task 변화는 성장보다 단조성 가정이 특히 불확실하다. [HSA 정의](https://academic.oup.com/nar/article/48/W1/W488/5815821); [실험 설계](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:107)
- 'w를 올리면 신호가 사라짐'은 선택한 요약이다. 추정량과 대비 표준오차가 함께 바뀌므로 nominal/BH 결과가 비단조이며, 한 가중치에서 발견된 q<.05를 검증된 용량 보정 결과로 해석할 수 없다. [코드](D:/Cander/School/3-2/system_bio/mtea/dose_weight_scan.py:33); [수치표](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:118)

### Gaps
- digitize 원좌표·반복별 데이터·불확실성 전파가 없어 사용자 성장 억제율 78%/92%, w=.7–.8의 정밀도는 검증하지 않았다. S1 Text는 대표 생물학적 반복임을 명시한다. [S1 816줄](D:/Cander/School/3-2/system_bio/mtea/si/flobak_S1.txt:816)
- Flobak main Figure 5는 SD error bar를 보여주지만, 이번 확인에서 biological replicate의 정확한 전체 n과 CI의 신뢰구간은 찾지 못했다. '독립 다중반복으로 통계적으로 확립' 같은 추가 문구는 피한다. [Figure 5](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004426)

## Q4. 24 h 전사체가 48 h 성장을 대리하지 못한다는 결론은 정당한가?

### Takeaway
'이 비교에서 증식 서명의 조건별 순위가 선행 성장실험과 일치하지 않았다'는 가능하다. '전사체 task 점수는 약효를 대리하지 못함'은 검증 범위를 넘는다. 권장: '현재 분석으로 task 점수와 성장효과의 연결은 검증되지 않았다.'

### Cited Findings
- `prolif_check.py`가 계산하는 것은 Tirosh S/G2M 유전자 평균의 조건별 발현차다. metabolic task 점수와 성장 간 회귀·예측성 검정을 수행하지 않는다. [prolif_check 4–19줄](D:/Cander/School/3-2/system_bio/mtea/prolif_check.py:4)
- 성장 값은 기존 그림에서 읽은 단일 대표 반복의 숫자를 print로 고정 입력한다. 같은 실험의 개별 RNA 샘플과 짝지어진 phenotype이 아니다. [prolif_check 37–39줄](D:/Cander/School/3-2/system_bio/mtea/prolif_check.py:37); [extension_report 124줄](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:124)
- Flobak은 48 h 성장 endpoint; Benedicto는 24 h RNA이며, Flobak 임피던스 cell index는 직접 세포수·세포사멸 측정과 구분된다. [Flobak Methods](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004426); [Benedicto Methods](https://www.nature.com/articles/s41540-025-00586-y)
- PI3Ki 단독 성장 억제율은 MEK 병용 실험에서 52%, TAK 병용 실험에서 73%로 다르다. 서로 다른 성장 실험의 값을 한 순위로 합친다. [보고서 111–122줄](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:111)
- Tsirvouli는 시작 시 80–90% confluency를 보고하고, 24 h DMSO에서도 cell-cycle 변화가 contact inhibition 등에 영향받을 수 있다고 토의한다. [Tsirvouli Methods·cell-cycle discussion](https://www.biorxiv.org/content/10.1101/2024.03.12.584561v1.full)

### Inferences
- S/G2M 평균 변화는 phase fraction, 세포주기 정지 위치, 적응·스트레스 반응의 영향을 받을 수 있다. 차후 누적 성장이 감소하는 조건에서 24 h 일부 증식 유전자가 올라가는 것은 약효 대리 가능성 전체를 기각하는 증거가 아니다. 직접 matched 시점·용량·반복 성장값이 없다. [분석 범위](D:/Cander/School/3-2/system_bio/mtea/prolif_check.py:14); [실험 차이](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:124)
- 성장 synergy 원논문 근거의 강도와 특정 metabolic mechanism의 근거 강도를 분리해야 한다. 원논문 CI 보고는 성장 phenotype을 뒷받침하지만 사용자 task의 인과 기여를 검증하지 않는다. [Flobak Results](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004426); [현 분석](D:/Cander/School/3-2/system_bio/mtea/prolif_check.py:37)

### Gaps
- 동시 RNA/성장, 24/48 h 양쪽, matched half-dose single 조건이 없어 대리성의 실패 원인과 예측 성능을 평가하지 못한다. 대규모 재분석은 하지 않았다. [명시된 설계](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:124)

## Q5. Tsirvouli 전량 투여·기울기·독립 검증 주장은 맞나?

### Takeaway
전량 투여는 기울기가 아니라 Methods로 확인 가능하다. 별도 실험 정황은 충분하나, 같은 연구팀·AGS 모델의 후속 데이터이므로 외부 독립검증으로 표현하지 말 것. 추가 중요 오류: 입력은 단순 log2-normalized 값이 아니라 VSN 변환값이며 차이를 정확한 log2FC로 부른다.

### Cited Findings
- bioRxiv Methods는 PI 0.70 µM, 5Z 0.50 µM, PD 0.035 µM를 단일 또는 병용으로 넣는다고 명시한다. '전량'은 해당 단일 GI50 기준이며 2025의 병용 반량과 다르다. [Tsirvouli Methods, Transcriptomics sample preparation](https://www.biorxiv.org/content/10.1101/2024.03.12.584561v1.full); [Benedicto Table 2](https://www.nature.com/articles/s41540-025-00586-y)
- 공식 repo `tr_limma_per_Timepoint.R`는 gene_counts_flt_norm_wCTRL.csv를 vsnMatrix/predict로 만들고, 별도 공식 differential expression은 TMM→voom→lmFit→eBayes로 수행한다. 사용자가 내려받은 counts_norm.csv는 전자 파일이다. [공식 코드](https://github.com/Eirinits/AGS_synergies/blob/main/Transcriptomics/scr/tr_limma_per_Timepoint.R); [사용자 코드 8, 24–42줄](D:/Cander/School/3-2/system_bio/mtea/tsirvouli_analysis.py:8)
- VSN은 큰 값에서 log2에 접근하는 분산 안정화 변환이며, 일반 log2와 동일하지 않다. 현재 평균차는 'VSN 변환 발현차'라 부르는 편이 정확하다. [VSN 공식 vignette](https://bioconductor.org/packages/release/bioc/vignettes/vsn/inst/doc/A-vsn.html)
- 공식 metadata는 96행(6처리×5시간×3반복+6untreated); 실제 matrix header는 94샘플(처리89+untreated5). PD_1_Rep1과 untreated_0_Rep1이 제외되어 PD 1 h는 n=2, 다른 처리·시간은 n=3이다. 공식 코드가 Sample_ID X101·X63을 제외한다. [공식 metadata](https://github.com/Eirinits/AGS_synergies/blob/main/Transcriptomics/Data/sample_info_to_use.tsv); [원 metadata](https://github.com/Eirinits/AGS_synergies/blob/main/Transcriptomics/Data/Raw_data/sample_info.tsv); [공식 DE 코드](https://github.com/Eirinits/AGS_synergies/blob/main/Transcriptomics/scr/tr_limma_per_Timepoint.R); [로컬 metadata](D:/Cander/School/3-2/system_bio/mtea/tsirvouli/sample_info.tsv:1)
- Tsirvouli는 HiSeq4000 single-end 64 cycles; Benedicto는 HiSeq2000 paired 2×64 bp, 24샘플 n=4. Tsirvouli raw filename에 2022 sequencing-run folder와 별도 Sample_ID가 있으며 연구팀 Flobak이 양 논문에 겹친다. [Tsirvouli Methods·저자](https://www.biorxiv.org/content/10.1101/2024.03.12.584561v1.full); [Benedicto Methods](https://www.nature.com/articles/s41540-025-00586-y); [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285616); [repo raw metadata](https://github.com/Eirinits/AGS_synergies/blob/main/Transcriptomics/Data/Raw_data/sample_info.tsv)
- 기울기는 |single|>.5, |partner|<.25로 유전자를 고른 뒤 원점회귀를 한다. 코드 출력은 slope1=same dose, .5–.7=reduced dose로 단정한다. [사용자 코드 60–76줄](D:/Cander/School/3-2/system_bio/mtea/tsirvouli_analysis.py:60)
- per-time 모델은 condition cell-means만 쓰며 rep를 batch로 모델에 포함하지 않는다. 파일명 rep가 같은 날 배치인지 확인해야 한다. [사용자 코드 13, 24–29줄](D:/Cander/School/3-2/system_bio/mtea/tsirvouli_analysis.py:24)

### Inferences
- 기울기 변화는 용량 외에 상호작용, 동역학, 공유 DMSO 추정오차, noisy x 선택, regression dilution, VSN 척도 차이로 생긴다. slope1은 동일 용량의 필요충분조건이 아니다. 특히 partner 단독이 작다는 이유로 조합에서 partner 작용이 없다고 가정할 수 없다. [실제 선택·회귀식](D:/Cander/School/3-2/system_bio/mtea/tsirvouli_analysis.py:67)
- 권장 표현: '동일 연구팀의 별도 전량 투여 시계열에서 방향 일치 여부를 탐색했다.' '독립 dose-matched 검증'은 2025와 병용 용량이 다르며, VSN차와 원논문 logFC도 다르므로 과하다. [두 설계](https://www.biorxiv.org/content/10.1101/2024.03.12.584561v1.full); [2025 설계](https://www.nature.com/articles/s41540-025-00586-y)
- 가장 작은 후속 검증은 공개 공식 voom logFC 출력으로 핵심 유전자 방향·기울기 수치를 확인하고, 이미 명시된 투여량을 기준으로 분석을 라벨링하는 것. raw 재분석 없이도 VSN→logFC 명명 오류를 해결할 수 있다. [공식 출력 경로](https://github.com/Eirinits/AGS_synergies/blob/main/Transcriptomics/Output/20230904_DEGs_perContrast_all.csv); [공식 코드](https://github.com/Eirinits/AGS_synergies/blob/main/Transcriptomics/scr/tr_limma_per_Timepoint.R)

### Gaps
- bioRxiv v1의 raw transcriptomics accession은 XX placeholder. 따라서 두 논문 간 생물시료가 전혀 재사용되지 않았다는 identity-level 검증은 미완료. 플랫폼·sample-ID·용량·반복이 달라 별도 실험 해석은 타당하지만 외부 기관 검증은 아니다. [Tsirvouli Data and code availability](https://www.biorxiv.org/content/10.1101/2024.03.12.584561v1.full)
- 반복 라벨의 batch 의미, 분산 안정화 변환에 직접 선형모형을 적용한 방법의 수치 오차·검정 교정은 따로 확인하지 않았다. 분석 코드는 수정하지 않았다. [현 모델](D:/Cander/School/3-2/system_bio/mtea/tsirvouli_analysis.py:24)

