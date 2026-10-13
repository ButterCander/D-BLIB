# 생물학적 주장, 저자 재현과 BH 검수

## Q6. DepMap 필수성을 약 반응 해석에 어디까지 쓸 수 있는가?

### Takeaway
DepMap은 약물을 처리하지 않은 배양 조건에서 유전자 knockout이 성장에 미치는 영향이다. 후보의 낮은 기준선 의존성은 확인되지만, 약물 처리 중 적응·시너지 기전의 불필요성을 증명하지 않는다.

### Cited Findings
- Chronos 공식 문서는 gene effect를 성공적 knockout에 따른 상대 성장률 변화로 정의한다. 0은 변화 없음, 음수는 생존도 감소다. 보통 공통 필수유전자 중앙값 −1, 비필수유전자 중앙값 0으로 정규화한다. `−0.5`는 본 분석의 실용 cutoff다. 공식 문서는 음성대조 기반 hit calling도 제공한다. — [Chronos 공식 README](https://raw.githubusercontent.com/broadinstitute/chronos/master/README.md); [로컬 cutoff 및 분석 코드](D:/Cander/School/3-2/system_bio/mtea/external_check.py:14)
- 실제 저장된 AGS 행 ACH-000880을 다시 읽었다. ASNS −0.056825, PHGDH −0.007509, PSPH +0.003811, ODC1 −0.151842, SRM −0.096478, TXNRD1 −0.206861. RRM1 −2.277380, RRM2 −1.997960, DUT −2.203274, DHODH −1.543278, UMPS −1.651420. 보고서 반올림 수치는 맞는다. — [DepMap 저장 데이터](D:/Cander/School/3-2/system_bio/mtea/external/crispr_gastric_rows.tsv); [공식 릴리스 파일을 읽는 코드](D:/Cander/School/3-2/system_bio/mtea/depmap_stream.py:14)
- **CoQ 비필수라는 묶음 결론은 내부 모순이다.** COQ2 −0.529604, COQ3 −0.502147, COQ5 −0.513672로 모두 본 분석의 −0.5 cutoff를 통과한다. 경계 수준 의존성이며, 다른 강조 후보와 한꺼번에 비필수로 묶으면 안 된다. — [저장 데이터](D:/Cander/School/3-2/system_bio/mtea/external/crispr_gastric_rows.tsv); [보고서 CoQ 수치](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:162); [모순인 결론](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:168)
- AGS `Model_24Q4.csv`의 FormulationID는 **F-12K + 10% FBS**, OnboardedMedia는 MF-028-001이다. 원 논문 조건은 Ham’s F12 + 5% FCS다. '다른 배지'보다 이 구체 차이를 쓰는 편이 정확하다. 다만 OnboardedMedia가 해당 CRISPR screen의 실제 배지였는지는 이 파일만으로 확정하지 못했다. — [모델 메타데이터](D:/Cander/School/3-2/system_bio/mtea/external/Model_24Q4.csv); [원 논문 Cell culture](https://www.nature.com/articles/s41540-025-00586-y)
- CCLE 원 자료는 928개 세포주의 225개 대사체 프로파일이다. 약물 처리 전 세포주 비교값을 PI103 병용의 대사 변화로 해석할 근거는 없다. — [Li et al. 2019 원 논문](https://www.nature.com/articles/s41591-019-0404-8); [분석에서 AGS 백분위를 계산하는 코드](D:/Cander/School/3-2/system_bio/mtea/external_check.py:49)

### Inferences
- 허용 문장: **“이 후보들은 DepMap의 해당 배양 조건에서 강한 기준선 성장 의존성을 보이지 않았다.”** 약물에 의한 부분 발현 감소와 유전자 완전 knockout의 효과가 같다는 가정은 불필요하다.
- 금지 문장: “후보는 성장에 관여하지 않는다”, “약물 시너지 기전이 아니다.” 치료 조건에서 새로 생기는 의존성, 영양 공급, 대체 경로·isoenzyme 보상, knockout 시점과 치료 시점 차이가 남는다.
- **핵산 합성은 성장과 연관된 후보 축**이다. AGS의 강한 기준선 의존성과 일부 전사 하향이 함께 관찰된다. 이는 PI103–MEK 병용 성장 억제의 매개기전이라는 인과 증거가 아니다. 공통 필수성·세포주기 변화가 대안 설명이다.
- `task_dep`는 집합 내 **최소 Chronos**로 task를 순위화한다. 집합이 클수록 매우 음수인 유전자를 포함할 가능성이 크다. task 크기·공유유전자를 통제하지 않은 MW p<0.001은 우선순위 참고치다. — 근거 코드: [external_check.py](D:/Cander/School/3-2/system_bio/mtea/external_check.py:37)

### Gaps
- 해당 CRISPR screen의 실제 배지·serum lot·assay duration 및 drug-conditioned dependency는 확인하지 못했다.
- ASNS/ODC1/SRM/PHGDH 등의 AGS에서 약물 처리 중 loss-of-function, 영양물 보충 rescue, 동위원소 flux 실험은 이 패키지에 없다.
- 29.2% 대 12.7%는 보고서/코드에서 확인했으나 새로 대규모 계산하지 않았다. 이는 task-essential 모델 개념의 전반적 enrichment이며 해당 약물 시너지 검증은 아니다.

## Q7. 재현·BH와 최종 주장 강도는 맞는가?

### Takeaway
점수 방향의 재현과 유의성의 재현을 나눠야 한다. 59 대 2의 BH 격차를 해결하기 전 “견고한 PI3K 주 동인”을 강으로 둘 근거는 부족하다. 23/48은 시너지 검정 통과 수가 아니다.

### Cited Findings
- **원 논문은 TIDE와 TIDE-essential 모두 10,000 permutations를 사용했다.** 저자 1,000 대 재현 10,000의 차이로 현재 격차를 설명하면 틀린다. 논문은 비유의 유전자 logFC를 0으로 마스킹하고 min/absmax GPR을 적용한다고 명시한다. — [원 논문 Methods](https://www.nature.com/articles/s41540-025-00586-y); [저장된 원문](D:/Cander/School/3-2/system_bio/mtea/si/page.html:1046)
- 공식 MTEApy `calculate_pvalue`는 `min(P(random≤score), P(random≥score))`다. 따라서 양측 p는 `min(1,2p)`이며 raw p<0.025와 양측 p<0.05의 비교는 일관된다. **+1 보정 없이 p=0을 허용**한다. 유효한 permutation p로 보려면 Monte Carlo uncertainty와 원래 null의 타당성은 따로 점검해야 한다. — [공식 소스](https://raw.githubusercontent.com/bsc-life/mteapy/main/src/mteapy/utils.py); [로컬 utils.py](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/utils.py:221)
- MOESM4 전체 187 task를 작은 메모리 계산으로 다시 읽고 `min(1,2p)`에 표준 BH step-up을 적용했다. TAKi 31→0, MEKi 19→0, PI3Ki 78→59, PI3Ki-TAKi 56→3, PI3Ki-MEKi 64→27. 표의 숫자는 재확인됐다. TAKi/MEKi 시트의 상단 표제는 뒤바뀌어 있다. — [Supplementary Data 3 로컬](D:/Cander/School/3-2/system_bio/mtea/si/MOESM4.xlsx); [보고서](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:19)
- MOESM5의 PI3Ki-MEKi 48개 및 PI3Ki-TAKi 14개 task ID를 MOESM4의 병용 p/q에 연결했다. 23/48, 1/14, q<0.1에서는 36/48, 8/14가 재확인됐다. **이 q는 병용 대비 control의 task 점수 p에서 나온다.** 상호작용 대비 또는 병용–각 단일약 차이의 p가 아니다. — [Supplementary Data 4 로컬](D:/Cander/School/3-2/system_bio/mtea/si/MOESM5.xlsx); [보고서가 명시한 계산 대상](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:27); [BH p를 읽는 코드](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:77)
- 원 논문 strong synergy는 두 단일약의 **유의성 마스킹 후 MS=0**, 병용의 MS≠0이라는 판정이다. 따라서 “한쪽 유의, 다른 쪽 비유의”가 조건 차이의 유의성을 보장하지 않는 문제가 본질이다. weak은 병용 점수가 단일약 점수보다 더 극단적인 휴리스틱이다. — [원 논문 Eq.1–2](https://www.nature.com/articles/s41540-025-00586-y); [원문 저장본](D:/Cander/School/3-2/system_bio/mtea/si/page.html:1006)
- 재현은 gene symbol 중복의 첫 행 유지, 합계 count≥10 필터, PyDESeq2, `~batch+condition`, shrinkage, ALLGENES 여부에 따른 null 모집단 선택을 사용한다. 저자는 R DESeq2 1.40.0·apeglm을 명시했다. 현재 저장본에 저자 exact gene-level DEA table을 대조한 결과는 없다. — [재현 코드](D:/Cander/School/3-2/system_bio/mtea/repro.py:21); [null 모집단 선택](D:/Cander/School/3-2/system_bio/mtea/repro.py:49); [원 논문 DEA Methods](https://www.nature.com/articles/s41540-025-00586-y)
- 같은 10,000 permutations 로그에서 전체 유전자 null의 raw 유의 task 수는 PI103/PD/5Z/PI103+5Z/PI103+PD = 44/16/22/29/34다. 모델유전자 null은 37/29/14/40/48로 바뀐다. null 모집단 민감성이 실제 존재한다. — [p10k.log](D:/Cander/School/3-2/system_bio/mtea/p10k.log); [TAG와 ALL 플래그](D:/Cander/School/3-2/system_bio/mtea/repro.py:17)
- PI103은 PI3K 외 mTOR/TORC1, DNA-PK 억제 활성도 가진다. 2007 primary pharmacology 논문에서 TORC1 및 DNA-PK 억제가 확인된다. 따라서 PI103에 반응했다는 결과를 순수 PI3K 경로 인과성으로 바꿀 수 없다. — [Raynaud et al. 2007 원 논문](https://aacrjournals.org/cancerres/article/67/12/5840/533102/Pharmacologic-Characterization-of-a-Potent); [별도 primary PI103/mTOR 연구](https://pubmed.ncbi.nlm.nih.gov/18548104/)

### Inferences
- **상한선 철회:** “낙관적 p에 기초한 참고 수치이며 FDR 5% 제어가 확립되지 않았다”로 바꿔라. 상관 보정으로 p가 전부 커지고 task가 일부만 남는다는 수학적 보장은 없다. 특히 TIDE 대 gene-union CAMERA-like는 통계량과 검정 모집단도 다르다.
- **23/48 재명명:** “저자 synergy-call 중 23개는 병용 대 control의 차이가 BH를 통과했다.” “23개 synergy가 보정 후 살아남았다”는 부정확하다.
- **strong 8개 재명명:** “그 8개 task의 병용 대 control 차이는 BH를 통과하지 않았다.” “strong synergy 상호작용 p가 탈락했다”라고 쓰면 다른 가설로 바뀐다.
- **재현 판정:** “task 점수의 방향·순위는 대체로 유사했으나, p/q와 유의 task 수는 재현되지 않았다.” 높은 score r이 p 보정을 인증하지 않는다. 순위 상관은 배경 모집단·유전자 마스킹·분산 또는 p 계산 차이를 놓친다.

| 기존 주장 | 권장 강도·문장 |
|---|---|
| 성장 시너지 실재, 강 | 기존 실험에서 성장 수준 synergy가 보고됨. 디지타이즈한 1대표곡선만으로 새 정량 판정을 확정하지 않는다. 성장 검수 담당의 원 실험 반복·CI 검증 결과를 우선 적용. |
| PI3K 억제가 견고한 주 동인, 강 | **중:** PI103 포함 조건의 공통 전사 반응이 뚜렷하다. task 유의성 재현과 PI3K 선택적 인과성은 미확인. |
| task 시너지 이 설계로 확립 불가, 중 | **타당:** 반수 단일약 arm 부재와 분석·검정력 한계에서 확립되지 않았음. 전역적으로 효과 없음이라는 결론은 불가. |
| task 점수는 성장 약효를 대리하지 못함, 중 | **약~중:** 현재 24 h signature가 기존 48 h 성장 순위를 재현하지 않았다. task 점수 일반의 surrogate 실패를 증명하지 않는다. |
| 강조 후보는 AGS 성장 필수 아님, 중 | **측정 범위에서는 중:** 지정한 후보의 DepMap 기준선 강한 의존성 없음. CoQ 예외, 약물 조건에서 비필수는 미확인. |
| 성장 관련 축 핵산 합성, 중~약 | **약:** 핵산 합성은 강한 기준선 의존성 및 일부 전사 감소를 보인 탐색 후보 축. 시너지 매개 인과성 미확인. |
| mRNA는 한 층, 중 | 생물학적으로 타당. 본 결과는 mRNA와 관측 phosphosite 변화의 낮은 동조를 보인 보조 근거. flux 변화 또는 독립성 증명은 아님. |

### Gaps
- 저자 입력 DEA, 모델 버전/부속 행렬, 유전자 universe, masking, shrinkage를 한 단계씩 고정해 59 대 2를 분해하는 감사가 아직 없다. 이번 검수에서 큰 재계산은 하지 않았다.
- 모든 task p가 보정 후 증가한다는 전제는 검증되지 않았다. “상한선” 사용 불가.

## Q8. 놓친 한계와 우선순위

### Takeaway
새 분석을 더하기 전에 재현 격차와 서로 다른 검정 가설을 정리해야 한다. 이후 기능 해석은 task 구조·CRISPR·phosphosite·약물 표적의 조건 차이를 좁히는 쪽이 효율적이다.

### Cited Findings
- TIDE 점수는 pFBA 참고해를 통해 고른 반응에 min/absmax로 투영한 logFC 평균이다. 현재 `taskgenes.py`는 그 반응에 등장하는 모든 유전자의 합집합을 만든다. 따라서 gene-union 평균검정은 GPR, 반응 반복 가중, AND/OR 비선형을 유지한 TIDE 점수 검정이 아니다. — [taskgenes.py](D:/Cander/School/3-2/system_bio/mtea/taskgenes.py:7); [공식 TIDE score 소스](https://raw.githubusercontent.com/bsc-life/mteapy/main/src/mteapy/tide.py); [로컬 TIDE score](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/tide.py:196)
- 실제 `task_genes.json` 187개는 고유 유전자 집합 181개다. X19/X20 등 동일집합 task가 존재하며, X153/X155/X157은 같은 64개 유전자로 구성된다. 단순 task 개수는 독립 생물학적 발견 수가 아니다. — [task_genes.json](D:/Cander/School/3-2/system_bio/mtea/out/task_genes.json)
- 원 논문 자체는 187 task 중 149 task에서 대안 최적 flux solution을 보고했다. TIDE-essential은 task 수행에 필요한 유전자이며 AGS의 실제 성장필수유전자 정의는 아니다. — [원 논문 TIDE-essential 및 FVA 결과](https://www.nature.com/articles/s41540-025-00586-y); [task-essential 행렬 정의](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/tide.py:23)
- phospho 분석은 같은 유전자의 site logFC를 단순 평균한다. 30분 phospho는 1 h RNA와 연결하며, “37% mRNA 변화 없음”은 PIPD 8 h에서 최대 |site logFC| site를 골라 **24 h mRNA |logFC|<0.25**로 정의한다. 동시점 검정이나 mRNA 변화 없음의 통계 검정이 아니다. — [시간 연결](D:/Cander/School/3-2/system_bio/mtea/phospho_analysis.py:22); [site 평균](D:/Cander/School/3-2/system_bio/mtea/phospho_analysis.py:47); [37% 계산](D:/Cander/School/3-2/system_bio/mtea/phospho_analysis.py:52)
- 효소 phosphosite enrichment는 site 단위 Fisher test다. 한 단백질의 여러 site가 중복관측으로 들어간다. `OR≈1`은 농축 증거 부재이지 대사효소 phospho 변화 부재가 아니다. — [phospho_analysis.py](D:/Cander/School/3-2/system_bio/mtea/phospho_analysis.py:32)
- 37% 요약·site 평균과 task 점수에는 단백질 abundance 보정 및 개별 site의 활성/억제 기능 annotation이 없다. 특정 site의 phosphorylation 증가가 효소 활성 증가를 뜻한다고 판단할 근거가 이 코드에는 없다. — [전체 phospho 분석](D:/Cander/School/3-2/system_bio/mtea/phospho_analysis.py:5); [원 보고서 한계](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:152)

### Inferences
- **최우선:** 저자 대비 재현 p 격차를 DEA/마스킹/모델/null 배경/점수/p 단계별로 분해. 발표 Q&A에서는 “점수 유사, 유의성 미재현”을 명시.
- **다음:** 저자 call, 병용 대 control q, 병용–단일약 차이, dose interaction을 서로 다른 질문으로 표기. BH 숫자를 synergy 입증 숫자로 재사용하지 않는다.
- **다음:** task의 고유 gene/reaction 구성과 중복·공유도를 요약하고, 최소 Chronos로 task 순위화한 결과는 크기 교란을 밝힌다.
- **후속 실험 후보:** 반수 단일약 arm과 동시점 성장/발현; 특정 후보의 약물 조건 CRISPR·rescue; 중요 phosphosite 기능 annotation 및 total-protein 보정; 실제 uptake/secretion 또는 동위원소 flux. 시간 없으면 새 분석 대신 이 순서를 Q&A 답변으로 남긴다.
- phospho 낮은 rho의 권장 표현: **“관측한 phosphosite 변화와 mRNA 변화는 낮은 단조 상관을 보였다.”** 독립층, 무관, mRNA가 기능을 반영하지 못한다는 표현은 강하다.
- PI103 표적 선택성, AGS 한 세포주 일반화, 배지 조건에서의 task feasibility/영양물 보충은 결과 범위를 좁힌다. 전사 감소가 대사 flux 감소 또는 성장 억제의 원인임을 자동으로 뜻하지 않는다.

### Gaps
- 모든 task의 flux해·essentiality를 원 배지·AGS expression에 맞춰 다시 구축하지 않았다.
- phospho primary dataset 전체 재분석, protein-normalized site 데이터 또는 phosphosite 기능 데이터베이스 조인은 하지 않았다.
- 이번 수치 확인은 저장된 출력·소규모 표/행 조회에 한정했다. 분석 원본을 편집하거나 대규모 재계산하지 않았다.
