# 통계·구현 독립 검수

범위: `camera_like.py`, `hsa_analysis.py`, `validate_checks.py`, `validate_interaction.py`, `nullde.py`, `inter.py`, `taskgenes.py`, 연결된 MTEApy 코드·출력. 원본 수정 없음. AGENTS.md 파일은 작업 경로 재귀 검색에서 없었고, 사용자 제공 지침 적용. ponytail full/review 읽음. 아래 '확인'은 코드·파일에서 확인한 사실, '추론'은 검수자의 해석. 원 분석 대규모 재실행 없음.

## Q1. CAMERA-like는 무엇을 검정하며, 공식 CAMERA·FRY·ROAST와 맞는가?

### Takeaway

잔차 평균 상관과 VIF를 직접 추정하는 방향은 타당함. 그러나 현재 구현은 공식 CAMERA의 검정과 다르며, TIDE의 p만 보정한 분석도 아님. '상관을 고려한 탐색적 gene-set 민감도 분석'으로 제시 가능. '보정된 확정 검정'으로 부르기는 어려움.

### Cited Findings

- **확인: 경쟁검정.** 로컬 검정의 분자는 `mean(z_in)−mean(z_out)`이다. 따라서 집합의 약물 반응 자체가 0인지가 아니라, 나머지 유전자보다 상대적으로 다른지를 검정한다. 전 유전자가 함께 변하면 해당 task가 실제 변해도 경쟁검정은 비유의일 수 있다. — [camera_like.py:34](D:/Cander/School/3-2/system_bio/mtea/camera_like.py:34), [CAMERA 원 논문, competitive null](https://pmc.ncbi.nlm.nih.gov/articles/PMC3458527/)
- **확인: rho 추정 식 자체는 CAMERA 원리와 맞음.** `max(0, 평균 pairwise residual correlation)`과 `VIF=1+(m−1)rho`를 쓴다. 평균을 제거한 OLS 잔차를 단위 길이로 정규화하는 방식은 공식 `interGeneCorrelation`의 QR 잔차 계산과 같은 원리를 따른다. — [camera_like.py:16](D:/Cander/School/3-2/system_bio/mtea/camera_like.py:16), [camera_like.py:24](D:/Cander/School/3-2/system_bio/mtea/camera_like.py:24), [공식 CAMERA 코드 L2–18](https://raw.githubusercontent.com/bioconductor-source/limma/master/R/geneset-camera.R)
- **확인: 자유도 처리 불일치.** 공식 CAMERA는 고정 rho일 때 gene 수 `G−2`, rho를 직접 추정할 때 `min(sample residual df,G−2)`를 사용한다. 이 설계의 residual df는 `24−9=15`. 로컬은 rho를 추정했어도 `G−2` 또는 `G_out−2`를 쓰므로 rho 추정 불확실성을 고정 rho처럼 취급한다. `G−2`가 언제나 잘못된 것은 아님. **추정 rho와 함께 사용하는 것이 문제**다. — [camera_like.py:38](D:/Cander/School/3-2/system_bio/mtea/camera_like.py:38), [validate_interaction.py:35](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:35), [validate_checks.py:54](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:54), [공식 CAMERA 코드 L71–79](https://raw.githubusercontent.com/bioconductor-source/limma/master/R/geneset-camera.R)
- **확인: 나머지 구현도 다름.** 공식 CAMERA는 moderated gene t→z, 집합과 배경의 pooled within-group variance를 쓴다. 로컬은 DESeq2 Wald z에 `var(z_out)`만 사용하며, rho는 별도 무가중 `log2(normalized count+1)` OLS 잔차에서 얻는다. 이는 표준 CAMERA와 동일한 분포 보장이 있는 조합이 아니다. — [camera_like.py:13](D:/Cander/School/3-2/system_bio/mtea/camera_like.py:13), [camera_like.py:36](D:/Cander/School/3-2/system_bio/mtea/camera_like.py:36), [공식 CAMERA 코드 L144–154, L193–198](https://raw.githubusercontent.com/bioconductor-source/limma/master/R/geneset-camera.R)
- **확인: rho=0.01은 진실값 아님.** 공식 문서는 fixed 0.01을 해석 가능한 순위와 엄격한 오류율 사이 절충으로 설명하고, `inter.gene.cor=NA`로 집합별 상관을 추정할 수 있다고 설명한다. 직접 추정은 원래 지원되는 방식이며, 0.01보다 크다는 사실만으로 추정치가 잘못됐다고 볼 수 없다. — [limma 공식 manual, camera p.45–46](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf)
- **확인: TIDE와 통계량·구조가 바뀜.** `taskgenes.py`는 task 반응에 연결된 모든 유전자의 중복 제거 union을 만든다. 원 TIDE는 반응별 GPR의 AND=min, OR=absmax를 계산한 뒤 반응 평균을 낸다. gene union은 복합체·동위효소 논리, 여러 반응에 반복 기여하는 유전자의 가중, 유전자 없는 반응의 0 기여를 버린다. — [taskgenes.py:7](D:/Cander/School/3-2/system_bio/mtea/taskgenes.py:7), [tide.py:196](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/tide.py:196), [utils.py:167](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/utils.py:167)
- **확인: 187→120은 크기 필터.** 최종 task은 측정 유전자 수 **5≤m≤200**. essential 검증 집합만 3≤m≤200. 읽기 전용 표준 라이브러리 점검: 원 task 187, 측정 유전자 16,932, 작은 task 56개·큰 task 11개 제외하여 120개. — [validate_interaction.py:18](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:18), [validate_checks.py:18](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:18), [검증 로그:1](D:/Cander/School/3-2/system_bio/mtea/validate_checks.log:1), [task_genes.json](D:/Cander/School/3-2/system_bio/mtea/out/task_genes.json)
- **확인: 동일 gene set 중복.** 120 task명 중 116개 유일 gene set. 동일 그룹은 X52/X139, X153/X155/X157, X192/X193. 이 사실만으로 BH가 틀리는 것은 아님. 다만 여러 task명이 독립된 생물학적 재현처럼 보이지 않도록 밝혀야 함. — [task_genes.json](D:/Cander/School/3-2/system_bio/mtea/out/task_genes.json), [validate_interaction.py:18](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:18)
- **확인: 배경 변경.** 원 재현에서 `ALLGENES`에 따라 metabolic-model 유전자만 쓰거나 전 발현 유전자를 쓴다. CAMERA-like는 전 발현 유전자가 배경. 배경은 경쟁검정의 질문을 바꾸므로 전 유전자와 발현 metabolic 유전자 배경을 같은 재현으로 취급하면 안 됨. — [repro.py:17](D:/Cander/School/3-2/system_bio/mtea/repro.py:17), [repro.py:51](D:/Cander/School/3-2/system_bio/mtea/repro.py:51), [camera_like.py:42](D:/Cander/School/3-2/system_bio/mtea/camera_like.py:42)
- **확인: FRY/ROAST는 self-contained.** 집합 내 변화 자체를 검정하고 회전으로 유전자 의존성을 보존한다. FRY는 ROAST의 빠른 근사. CAMERA와 같은 귀무가설을 검정하는 대체품이라고 표현하면 안 됨. — [limma 공식 manual, roast p.219–221](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf)

### Inferences

- **우선순위는 적절함:** 다중검정→상호작용 대비→상관 민감도→음성·양성 대조. 다만 마지막 단계에서 통계량을 바꿨으므로 '같은 TIDE 상호작용 검정을 교정해 41→0'보다는 '다른 상관 인지형 gene-set 검정에서 0'이 정확함. — [inter.py:53](D:/Cander/School/3-2/system_bio/mtea/inter.py:53), [validate_interaction.py:26](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:26)
- **BH 0/min q=.075의 안전한 해석:** '현재 단순화한 120 task gene-union 경쟁검정에서 BH 5%를 통과한 집합은 없었다.' 이 값만으로 실제 metabolic task interaction이 없음, 공식 CAMERA에서도 없음, TIDE p가 교정됐음은 입증되지 않음. — [extension_report.md:98](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:98), [validate_interaction.py:35](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:35)
- **가장 작은 추가 검증:** 같은 사전 지정 gene universe·set 정의·설계·대비로 limma-voom `camera(inter.gene.cor=NA)`를 한번 비교. 집합 자체의 비가법성 질문에는 FRY/ROAST를 별도로 추가. 이 검정들 역시 용량 설계 누락이나 mRNA→flux 해석을 해결하지는 못함. — [공식 CAMERA 코드](https://raw.githubusercontent.com/bioconductor-source/limma/master/R/geneset-camera.R), [limma 공식 manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf)
- **저자 BH 카운트는 엄밀한 상한선 아님.** 상관 무시로 낙관적일 수 있다는 방향은 맞지만, 통계량·배경·유전자 구조가 달라 '교정하면 반드시 이보다 적다'는 수학적 단조성은 없음. '상관을 고려하지 않은 p에 BH를 적용한 참고 결과' 권장. — [extension_report.md:103](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:103), [CAMERA 원 논문](https://pmc.ncbi.nlm.nih.gov/articles/PMC3458527/)

### Gaps

- 5–200 필터를 결과를 보기 전에 고정했는지는 기록 없음. 크기 기반 필터 자체를 사후 p 선택으로 단정할 근거도 없음. 제외된 67 task를 결론 범위에서 명시할 필요. — [validate_interaction.py:19](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:19)
- 공식 CAMERA/FRY/ROAST 직접 실행 없음. 본 검수에서 df15·pooled variance 소규모 수치 비교를 시도했으나 사용한 bundled Python에 scipy가 없어 중단. 설치·대규모 재분석 안 함. 현재 최소 q=.075는 문서 기록을 확인했고 독립 재계산은 안 했음. — [extension_report.md:98](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:98)

## Q2. max p HSA, 방향 선택, 공분산, TOST

### Takeaway

사전 고정한 방향에서 '병용이 두 단일약을 모두 넘어서는가'는 IUT의 좋은 용도이며 `max(one-sided p)`가 맞음. 현재 구현은 signed one-sided p를 보존하지 않고 양측 p·사후 방향 필터를 쓰며, 원 TIDE 점수 차 자체도 검정하지 않음. 'HSA 향상 미검출'은 유지 가능, '동등'은 불가.

### Cited Findings

- **확인: 실제 p는 max 양측 p.** MTEApy p는 `min(lower-tail,upper-tail)`이다. `hsa_analysis.py`는 `min(1,2*max(pA,pB))`를 계산하므로 max(two-sided p)와 같고, 점수 sign으로 up/down/mixed를 구분한다. mixed도 BH에 넣은 뒤 마지막에 제외한다. — [utils.py:238](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/utils.py:238), [hsa_analysis.py:17](D:/Cander/School/3-2/system_bio/mtea/hsa_analysis.py:17), [hsa_analysis.py:20](D:/Cander/School/3-2/system_bio/mtea/hsa_analysis.py:20)
- **확인: 유전자 대비 계산은 올바름.** DESeq2에 계수 대비 벡터를 전달한다. 공유 control을 단순 독립 추정치처럼 빼는 방식이 아니며, 계수 대비 대수는 ~1e−15 정확도 확인. 다만 이 검증은 gene LFC의 대수만 보장한다. task score 대수·SE·HSA p의 정확성은 보장하지 않는다. — [inter.py:35](D:/Cander/School/3-2/system_bio/mtea/inter.py:35), [validate_checks.py:66](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:66), [검증 로그:4](D:/Cander/School/3-2/system_bio/mtea/validate_checks.log:4)
- **확인: 원 TIDE score HSA 아님.** 대비별 LFC를 DEG masking한 다음 GPR로 task score를 만든다. GPR AND=min/OR=absmax는 비선형이고 masking도 대비별로 바뀐다. 따라서 `TIDE(AB−A)`는 `TIDE(AB)−TIDE(A)`와 일반적으로 다르다. TIDE-essential도 masking 전에는 선형 평균이지만 대비별 masking 후에는 점수 차와 일치하지 않는다. — [inter.py:53](D:/Cander/School/3-2/system_bio/mtea/inter.py:53), [utils.py:112](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/utils.py:112), [utils.py:167](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/utils.py:167)
- **확인 가능한 반례:** 한 OR 반응에서 AB 유전자 값=(1,2), A=(2,1)이면 absmax(AB)=absmax(A)=2라 원 반응 점수 차는 0. 하지만 absmax(AB−A)=absmax(−1,1)=−1(현재 tie 방식)이 된다. 이는 대비 대수 검증으로 잡히지 않는다. — [utils.py:120](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/utils.py:120), [tide.py:196](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/tide.py:196)
- **확인: 동등성은 비유의 차이와 다름.** TOST는 생물학적으로 정한 동등성 구간 [−Δ,+Δ]를 먼저 정하고 두 단측 검정을 수행한다. α=.05일 때 대응하는 90% CI가 구간 안에 들어야 함. Δ는 검출 실패 후 유리하게 고르면 안 됨. — [Lakens 저자 설명, Equivalence Testing](https://lakens.github.io/statistical_inferences/09-equivalencetest.html), [Lakens 2017 원 논문](https://doi.org/10.1177/1948550617697177)

### Inferences

- **IUT 공분산:** 사전 방향 s에서 `δA=s(AB−A)`, `δB=s(AB−B)`라 하면 H0=`δA≤0 또는 δB≤0`, H1=`δA>0 그리고 δB>0`. `pHSA=max(pA,pB)`의 유효성은 두 p의 독립성을 요구하지 않음. 참인 구성요소 하나라도 있으면 전체 기각이 그 구성요소 기각을 요구하기 때문. 공유 AB로 인한 공분산은 각 contrast SE를 계산할 때와 검정력 평가에 필요. — [IUT 원 논문 Berger·Hsu 1996](https://doi.org/10.1214/ss/1032280304), [inter.py:38](D:/Cander/School/3-2/system_bio/mtea/inter.py:38)
- **현재 max(two-sided)는 무조건 틀렸다고 단정하면 안 됨.** 대칭·유효한 p와 올바른 방향 아래에서는 fixed-direction one-sided IUT보다 보수적일 수 있음. 그러나 현재 코드는 선택된 tail 방향을 버렸고, sign 기준은 raw score=0인데 competitive null은 배경 값에 의해 중심이 달라질 수 있다. 이 조건들이 입증되지 않아 '정식 directional HSA'로 인증하기 어렵다. — [hsa_analysis.py:17](D:/Cander/School/3-2/system_bio/mtea/hsa_analysis.py:17), [utils.py:238](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/utils.py:238)
- **두 방향을 탐색할 때의 단순 수정 원리:** 각 방향의 signed one-sided p로 IUT를 만들고 방향 둘의 union을 보정: `min(1,2*min(max(pA_up,pB_up),max(pA_down,pB_down)))`. 대칭 null·같은 sign 조건이면 max(two-sided)와 대응할 수 있지만 mixed는 유효한 HSA p=1로 처리해 BH 전에 넣어야 한다. 반대 방향의 강한 변화도 작은 `pIU`가 되게 두고 사후 제거하면 그 p가 다른 task의 BH threshold를 바꾸므로 해당 HSA FDR 보장을 따로 증명해야 한다. — [hsa_analysis.py:18](D:/Cander/School/3-2/system_bio/mtea/hsa_analysis.py:18), [IUT 원 논문](https://doi.org/10.1214/ss/1032280304)
- **최소 수정안:** HSA 방향·효과 척도·생물학적 유익성을 먼저 정하고, 같은 척도의 AB−A/AB−B 단측 검정→max p→task 전체 BH. 'more extreme'가 대사적으로 더 유익한지는 별도 가정이며 up/down 모두 자동 synergy로 부를 수 없다. 반수 용량 병용을 전량 single과 비교하는 현재 HSA는 더 엄격한 benchmark지만, 반수 용량 single에 대한 표준 HSA를 대체하지 못한다. — [inter.py:16](D:/Cander/School/3-2/system_bio/mtea/inter.py:16), [hsa_analysis.py:18](D:/Cander/School/3-2/system_bio/mtea/hsa_analysis.py:18)
- **TOST 현재 적용 불가:** task별 효과 크기 및 sample-based CI가 없음. competitive mean-z−background 검정의 비유의는 원 TIDE 단위 동등성을 의미하지 않는다. 뒤의 0.2–0.6 '해상력'도 Δ 근거로 쓰면 안 됨. — [validate_interaction.py:30](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:30), [validate_checks.py:125](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:125)

### Gaps

- 원 reaction-level TIDE score의 AB−single 차·CI·directional HSA는 미검정. 현재 저장된 tail-min p만으로 특정 방향의 단측 p를 항상 복원할 수 없음. — [hsa_analysis.py:9](D:/Cander/School/3-2/system_bio/mtea/hsa_analysis.py:9)
- 적절한 동등성 한계 Δ를 결정할 외부 생물학적 검증·사전 계획 없음. 따라서 '동등'보다 '추가 효과를 검출하지 못함'이 안전. — [extension_report.md:101](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:101)

## Q6. null·spike-in은 실제 파이프라인을 검증했는가? FPR·검정력·BH 주장

### Takeaway

자기 점검은 좋은 탐색적 검증. 0.387 대 0.037, 0.3 log2 spike 검출 80%/BH42%는 로그에서 확인됨. 그러나 최종 DESeq2/TIDE/HSA의 오류율·검정력이 검증됐다고 전이할 수 없다. 0.2–0.6을 실험 설계의 검출 한계로 부르는 근거는 특히 부족함.

### Cited Findings

- **확인: 수치 기록은 존재.** 30개 batch 내 라벨 셔플×5 single-control 대비=150 대비. task 120집합에서 gene-resampling FPR=.387, CAMERA-like=.037. 0.3 log2 spike의 CAMERA-like raw 검출=.80, BH=.42. — [검증 로그:16](D:/Cander/School/3-2/system_bio/mtea/validate_checks.log:16), [검증 로그:19](D:/Cander/School/3-2/system_bio/mtea/validate_checks.log:19), [검증 로그:27](D:/Cander/School/3-2/system_bio/mtea/validate_checks.log:27)
- **확인: 로그 설명 'Wald z'는 코드와 다름.** null/spike에서 unmasked `log2(normalized counts+1)` OLS fit, residual df15 t를 CDF→normal quantile로 바꾼 z를 사용한다. 최종은 DESeq2 Wald stat. 둘의 상관·SE·분포가 같다는 검증은 없음. — [validate_checks.py:33](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:33), [validate_checks.py:40](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:40), [validate_interaction.py:27](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:27), [extension_report.md:95](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:95)
- **확인: 'as MTEApy'는 실제 MTEApy 검정 아님.** 검증은 unmasked z를 복원 추출한 gene-mean null을 중심화하고 absolute tail에 +1 보정을 쓴다. MTEApy는 masked LFC를 비복원 shuffle, reaction GPR로 투영, 최소 단측 tail/count로 계산한다. '일반 gene-independence null의 취약성'은 지지하지만 원 TIDE FPR=.39를 측정한 것은 아님. — [validate_checks.py:48](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:48), [validate_checks.py:55](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:55), [tide.py:253](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/tide.py:253), [utils.py:238](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/utils.py:238)
- **확인: full pipeline null은 별도 1 seed.** `nullde.py`는 실제 DESeq2→mask→TIDE지만, 한 셔플의 5 대비에서 DEG가 0이라 task score null이 퇴화했다. 실제 파이프라인 전체 오류율 측정의 충분한 반복이 아니다. — [nullde.py:6](D:/Cander/School/3-2/system_bio/mtea/nullde.py:6), [nullde.py:29](D:/Cander/School/3-2/system_bio/mtea/nullde.py:29), [extension_report.md:94](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:94)
- **확인: spike는 synthetic single-control shift.** 60 trials마다 같은 실험 데이터의 라벨을 다시 섞고 task 하나의 모든 유전자에 같은 shift를 한 조건 4샘플에 더함. 원 gene covariance·기저 전처리를 공유하는 내부 스트레스 테스트이며, 실제 반수 용량 interaction나 두 대비 HSA의 power 실험은 아니다. — [validate_checks.py:112](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:112)
- **확인: 해상력 0.2–0.6의 계산 대상이 다름.** `validate_checks.py`는 AB−PI3Ki 대비의 padj≥.05 LFC를 0으로 바꾼 뒤 random gene mean의 95% quantile을 구한다. 이 값은 실제 task별 biological-replicate sampling SE도, GPR task score CI도, 정해진 power의 minimum detectable effect도 아니다. — [validate_checks.py:126](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:126), [검증 로그:31](D:/Cander/School/3-2/system_bio/mtea/validate_checks.log:31)
- **확인: 원 permutation p는 0 허용.** `min(count_low/B,count_high/B)`여서 count=0이면 p=0. permutation p에 관측 permutation을 포함하지 않을 때 0 값·과소추정은 다중검정에서 특히 문제라는 원 논문이 있다. — [utils.py:238](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/utils.py:238), [Phipson·Smyth 2010 저자 원문](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf)

### Inferences

- **calibrated의 범위:** '이 데이터에서 생성한 OLS synthetic null에 대해 raw FPR=.037로 관찰됐다'가 정확함. p 값을 조정하거나 독립 holdout에서 검증한 것은 아님. 세밀한 정규화·필터 선택을 포함한 full pipeline이 교정됐다는 표현은 철회 권장. — [validate_checks.py:90](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:90), [validate_interaction.py:1](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:1)
- **반복 독립성 주의:** 150 대비가 150개 독립 실험은 아님. 같은 셔플 내 5 대비는 공유 control·샘플로 의존하고 task도 겹침. FPR의 분모 18,000개를 독립 Bernoulli로 취급하면 CI가 부당하게 좁아질 수 있음. shuffle 단위 집계·불확실성 보고가 적합. — [validate_checks.py:92](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:92), [validate_checks.py:97](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:97)
- **해상력 문구 철회 권장:** '작은 효과에 대한 검정력이 제한될 가능성'은 설계·반복 수로 합리적이나 'task당 약0.2–0.6 log2 이하는 이 설계로 검출 불가'는 현재 계산으로 증명되지 않음. TOST Δ나 저자 0.08 마진의 검출 불가능 증거로 쓰지 말 것. — [extension_report.md:97](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:97), [validate_checks.py:129](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:129)
- **1000→10000은 MC 정밀도 개선.** '1000이면 BH가 지나치게 엄격'은 일반 법칙이 아니다. 187 tests/α=.05이면 rank1 cutoff=.000267, rank8=.002139. 현재 min-tail p의 비영점 최소 양측 p=.002는 rank≥8에서 통과 가능하며 p=0은 오히려 낙관적. 표준 +1 보정이면 최소 양측 p≈.002 역시 고립된 최상위 task를 막지만 충분히 많은 작은 p는 통과한다. — [extension_report.md:17](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:17), [utils.py:238](D:/Cander/School/3-2/system_bio/mtea/mteapy/src/mteapy/utils.py:238), [Phipson·Smyth 2010](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf)
- **BH 설명의 최소 조건:** 조건·method별 task family 안에서의 보정은 보고된 카운트의 정의로 합리적. 그러나 여러 조합·대비·시간·방법·용량 가중치를 뒤진 뒤 최선 결과를 고르면 전체 탐색의 FDR은 통제되지 않음. 겹치는 task는 그 자체로 BH 금지 근거가 아니며, task 수를 독립 경로 수로 읽지 않는 것이 핵심. — [validate_interaction.py:36](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:36), [extension_report.md:115](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:115)
- **최종 문구 권장:** '원 연구의 task 수준 synergy는 용량 설계와 통계적 검증의 한계로 확립됐다고 보기 어렵다. 추가 분석에서도 상관 인지형 gene-set 경쟁검정의 BH 유의 신호는 제한적이었다. 이는 synergy 부재·PI3Ki와 동등성의 증명이 아니다.' — [extension_report.md:98](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:98), [extension_report.md:101](D:/Cander/School/3-2/system_bio/mtea/extension_report.md:101)

### Gaps

- null/spike 각 trial의 원 p 배열이 저장되지 않아 FPR·CI·KS 독립 재계산 불가. 로그 수치와 코드 구조까지 확인. 평균 BH 기각수=.02는 FDR=.02라는 뜻도 아님. — [validate_checks.py:101](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:101), [validate_checks.py:107](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:107)
- 실제 DESeq2 Wald z 기반 final contrast의 null 오류율·HSA IUT·interaction power 검증 없음. masked original TIDE pipeline을 nondegenerate synthetic count data에서 검증하거나 공식 도구와 일치 확인이 필요하나 본 검수에서는 실행 안 함. — [validate_checks.py:40](D:/Cander/School/3-2/system_bio/mtea/validate_checks.py:40), [validate_interaction.py:27](D:/Cander/School/3-2/system_bio/mtea/validate_interaction.py:27)

