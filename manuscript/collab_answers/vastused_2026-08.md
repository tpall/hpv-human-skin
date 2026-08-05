# Vastused kolleegide kommentaaridele, august 2026

Kommentaarid saabusid pärast 2026-08-02 esitlust. Numeratsioon järgib esitatud
järjekorda. Kõik arvud on jooksvast analüüsist (kohort 10 849 teeki pärast
SRP572392 väljajätmist), mitte esitluse hetkeseisust.

---

## 1. E1 ja E2 peaksid tulema samast transkriptist; rakuliinidel E1 >> E2, normaalses nahas ja kondüloomides vastupidi

**Tähelepanek peab paika ja osutus kõige väärtuslikumaks kommentaariks.**
Kvantifitseeritud: `bin/e1e2_integration_marker.R` → `results_e1e2/`.

| tier | n | mediaan E1:E2 | teeke E1>E2 |
|---|---:|---:|---:|
| rakuliin | 49 | **129×** | 44/49 |
| valimata nahk | 55 | 0,84× | 24/55 |
| produktiivsed kontrollid (tüükad/kondüloomid) | 41 | **0,27×** | 5/41 |
| EV | 10 | 0,30× | 2/10 |
| peri-lesionaalne | 15 | 0,16× | 0/15 |
| neoplaasia | 18 | 0,04× | 1/18 |

Rakuliin vs produktiivne kontroll: **AUC 0,968**, Wilcoxon *P* = 2,8×10⁻¹⁴.

**Mehhanism on kinnitatud sisemise kontrolliga.** HaCaT31 (SRP560153, 6 teeki)
on rakuliin, mis hoiab HPV31 genoomi **episomaalselt** — tema suhe on 0,6–1,0×,
täpselt nagu primaarkoel, mitte 129×. Rakuliin ilma integratsioonita käitub nagu
kude, seega marker loeb E2 ORF-i terviklikkust, mitte kultuuris kasvamist.

**Praktiline väärtus:** see on metaandmetest sõltumatu saastemarker. Ta leidis
**28 teeki**, mis on integratsioonimustriga, aga rakuliiniks märkimata — nende
hulgas `A375 shNIPBL`, `SKMEL-5 shHOXD13`, `624Mel_siMITF` (melanoomiliinid,
mille `tissue_source` on „skin") ja UT-SCC cSCC liinid siRNA-ga (diagnoosiks
„cutaneous SCC"). Selle tulemusena lisati **`engineered` tier** (226 teeki:
159 nahast, 36 neoplaasiast, 31 mujalt), mis teravdas peakontraste
(M2b neoplaasia vs nahk 9,1 → 12,3).

**Aus piirang:** integreerunud rühm on peaaegu tervenisti HPV18/HeLa (34-st 32
rakuliini), seega tüüp ja integratsioon on osaliselt segunenud. HaCaT31 lahutab
„rakuliin = kõrge suhe" seose, aga mitte „HPV18 = kõrge suhe" oma. Lõplikuks
lahutamiseks oleks vaja integreerunud HPV16 liine (SiHa, CaSki); meil ei ole
ühtegi piisava kattega, et tüüpi määrata.

---

## 2. Kui saaks vähkide metaandmeid sorteeritult näha

Tehtud: **`results_comments/neoplasia_metadata_sorted.tsv`**, 406 rida, RPM
järgi kahanevalt. Veerud: SRR, uuring, pealkiri, koeallikas, diagnoos,
koekategooria, domineeriv HPV tüüp, lugemite arv, RPM, üle-läve lipp.

Kokkuvõte: 406 vähiteegist ületab tuvastusläve **10**, tüübimääratlus õnnestub
**9-l**. Tipus on SRR19520029 (HPV16, 765 RPM), SRR19520005 (HPV12, 523),
ERR1971044 ja ERR2274877 (mõlemad HPV182, ~500). Ehk kogu neoplaasia-tier'i
„signaal" tuleb kümnest teegist 406-st.

---

## 3. HPV nakatab ainult mitmekihilist epiteeli; suuõõne vähk on tõusuteel ja HPV+

Bioloogia poolel on see täpselt see eeldus, millel põhineb uus
**mitte-permissiivse koe kontroll** (334 teeki rasv-, vere- ja närvikude:
28,7 % jälgi vs naha 27,8 %, 0 % üle läve). Vt käsikirja tulemuste osa.

Suuõõne kohta: kohordis on **39 suuõõne/neelu teeki**, neist 18 rakuliinid
(TIGK igemekeratinotsüüdid, HPV16 — juba teadaolevad saasteuuringud) ja
**21 primaarkude**. Viimased jagunevad kaheks ja käituvad täpselt nagu nahk:

| uuring | n | üle läve | tüübid | kude |
|---|---:|---:|---|---|
| SRP516430 | 9 | **0** | — | normaalne suulimaskest |
| SRP418337 | 3 | **3** | HPV6, HPV11 | kõri papilloom |
| SRP495254 | 9 | 2 | HPV6 | — |

Normaalne limaskest vaikib, papilloom on produktiivne — sama muster kui
nahk vs soolatüügas, ainult teises nišis.

**Aga kohordis EI OLE orofarüngeaalset kartsinoomi.** Tõusev HPV+ OPSCC, millele
kommentaar viitab, ei ole esindatud: meil on normaalne limaskest ja healoomuline
papilloom, mitte vähk. See on konkreetne laienduskoht — täpselt sama kuju nagu
anaalne pilootuuring, mis andis meile puuduva „kartsinoom vs kartsinoom"
kontrolli.

---

## 4. Mis on täpselt peri-lesional mucosa?

15 teeki kahest seotud uuringust (SRP314958 RNA-seq + SRP314961 MeRIP),
proovid `CON_1`–`CON_5`: **kondüloomiga patsiendi kondüloomikõrvane, tervena
näiv limaskest**. Sama patsient, sama sekveneerimisjooks, kahjustuse kõrvalt.

Analüüsis on nad omaette tier, sest kõik 15 on üle tuvastusläve (100 %) —
mediaan 39 RPM. See on kohordi kõige veenvam „kohalolu ilma nähtava
kahjustuseta" tõend ja seetõttu ei tohi nad olla valimata naha nulltaustas.

---

## 5. Kui palju oli HPV+ proove üldse? (~25 % nahas)

Mälu on lähedal, aga oluline on kaks eri arvu:

| | kogu kohort | valimata nahk |
|---|---:|---:|
| vähemalt üks HPV-lugem | **29,1 %** (3 154/10 849) | **28,0 %** |
| üle tuvastusläve (2,6 RPM) | **2,2 %** (235/10 849) | **0,1 %** |

Vahe nende kahe rea vahel ongi kogu töö sisu: jälgi leidub peaaegu kolmandikus
teekidest, aga tuvastuslävele jõuab tuhandik. Ja mitte-permissiivne kude annab
sama 28,7 % jälgi, mis tõestab, et ülemine rida ei mõõda infektsiooni.

---

## 6. Beeta-HPV nakkab suhteliselt kiiresti, st on väga nakkav. Kus/kuna on aken?

Seda küsimust **ei saa selle andmestikuga vastata** ja see väärib ausat
piiritlemist. Meie kohort on läbilõikeline avalik RNA-seq: ei ole vanust,
ei ole korduvproove, ei ole seroloogiat. Omandamise aken (kirjanduses
varane lapsepõlv, perekonnasisene ülekanne, karvanääpsud reservuaarina)
nõuab longitudinaalset või vanuse-stratifitseeritud disaini.

Mida meie andmed sellele siiski lisavad: β/γ reservuaar on
**transkriptsiooniliselt vaikiv igas eas ja igas nišis, kus me teda näeme**,
seega „nakkavus" ja „aktiivsus" on eraldi teljed. Kui aken on varane, siis
see, mida meie mõõdame, on aknajärgne püsiseisund.

---

## 7. Kas oleks võimalik teha ka proteoomi analüüs?

Põhimõtteliselt jah, praktiliselt on kaks tõsist takistust.

**Andmete kättesaadavus.** Avalikke inimese proteoome (PRIDE, MassIVE) on
suurusjärke vähem kui RNA-seq teeke, ja nahaproove nende hulgas veel vähem.
Meie jõud tuleb 12 200 teegist; proteoomi poolel jääks n kolmekohaliseks.

**Tundlikkus.** Meie põhitulemus on *madala taseme jälg* — mediaan alla
0,01 RPM. Massispektromeetria avastamislävi viirusvalgule kudedes on selle
kohta suurusjärke kõrgem, seega negatiivne proteoomitulemus ei eristaks
„puudub" ja „allpool läve" vahel. Proteoom oleks tugev **positiivsetes**
kontrollides (soolatüügas, EV — kus L1 on rikkalik) ja seal kinnitaks
produktiivsuse valgu tasemel. Nahavähis, kus meie väide on „vaikiv",
proteoom väidet ei lahendaks.

Soovitus: kui teha, siis suunatult — L1/E7 valgud produktiivsetes
kontrollides vs cSCC, mitte avastuslik proteoomika.

---

## 8. Kui genoom on integreerunud, siis on E2 ORF puudu; kuidas see DNA andmebaasist paistaks?

Loogika on õige ja see on täpselt sama mehhanism, mis annab küsimuse 1 tulemuse.
DNA poolel oleks allkiri: kontiig, millel L1 ja E6/E7 on olemas, aga **E2 ORF
puudub või on katkenud**, pluss lugemid integratsioonipiiril (viirus–peremees
kimäärsed lõigud).

**Shen'i andmestikuga seda testida ei saa, ja põhjus on struktuurne.** Nende
kataloog on üles ehitatud **L1-le**: HMM-otsing, klasterdus ja tüpiseerimine
käivad kõik L1 järgi (`feb_7_L1_nuc_hmmer_env.tsv`, `pv_all_L1.domtbl`,
`all_full_l1s_logan_and_pr_90.uc`). Avaldatud väljundites ei ole kontiigide
järjestusi ega ORF-tasemel annotatsiooni, ainult klastri- ja taksonoomiatabelid.
Seega on see andmestik E2 staatuse suhtes konstruktsiooni poolest pime.

Testimiseks oleks vaja Logani kontiige endid ja ORF-annotatsiooni; PaVE GFF
kaheksa ORF-iga on meil juba olemas (`hpv_gene_annotations.gff`). See on
teostatav, aga eraldi töö.

Kõrvalmärkus, mis on oluline ka DNA-reservuaari võrdluse tõlgendamisel:
kuna Shen'i kataloog on L1-põhine, kirjeldab „kohalolu" seal L1 kohalolu,
mitte tervet genoomi.

---

## 9. Produktiivne HPV: vähid vs mitte-vähid anaalis

164 produktiivsest anaalteegist:

| | n |
|---|---:|
| kasvaja / SCC / CIN | **115** |
| kontroll | **1** |
| määramata (pärinevad anaal-SCC uuringutest) | 48 |

Ehk **mitte-vähi võrdlusrühma sisuliselt ei ole**. See on juba esitluse
piirangureal ja on põhjus, miks emakakaela poolel tellisime eraldi
469 kihistatud normaalset/kontrollteeki (töö 500834) — just selle lähtetaseme
saamiseks, mis anaalpiloodil puudus.
