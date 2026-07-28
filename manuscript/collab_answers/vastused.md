# Vastused kaasautorite küsimustele

Kõik arvud on arvutatud käesolevast andmestikust skriptidega `make_tables.R` ja `figure_genus.R`
(samas kaustas). Tekkinud tabelid: `tableA_assigned_type_calls.tsv`, `tableB_type_reads_orf.tsv`,
`tableC_coinfections.tsv`, `tableD_hpv_positive_cell_lines.tsv`; joonis: `figure_genus_prevalence.png`.

Kolme küsimuse vastus toob välja probleeme, mis on olulisemad kui küsimus ise. Need on koondatud
lõppu ("Leiud, mis nõuavad käsikirjas parandust").

---

## 1. Mis asjad on "engineered libraries"?

Need on avalikud RNA-seq raamatukogud, mille metaandmed näitavad, et rakke on enne sekveneerimist
geneetiliselt manipuleeritud. Klassifikaator on `bin/cell_line_patterns.py`, muster-loend
`ENGINEERED_PATTERNS`. Tabatakse:

- shRNA / siRNA / sgRNA, sh- ja si-prefiksiga sihtmärginimed (shTOP2A, siSASH1, shNC)
- knockdown / knockout, CRISPR, Cas9 / dCas9
- transduktsioon, transfektsioon, lentiviirus, stabiilsed transfektandid
- plasmiidiselgroog vektorites (pLVX, pLKO, pBABE, pCDH, pCMV jne)

Näiteid tegelikest kirjetest: `N/TERT-2G keratinocytes, APOE knockdown`, `SKMEL-5 shHOXD13-3`,
`BJ sgRNA-PURPL#1`, `melanocyte, sphere assay, siSASH1`, `RNAseq_501mel shNIPBL_#1`.

Erapooletus sõelal (11 556 raamatukogu) oli 198 raamatukogu "engineered" (need 12, mis olid
korraga ka rakuliin, loeti rakuliiniks, sest rakuliini lipp on ülimuslik). HPV+ oli neist 55,
positiivsuse määr 11,1% (95% UI 7,45 kuni 16,2) vs kliinilise koe 2,15%. Domineerivad tüübid:
HPV18 (22 raamatukogu), HPV31 (4), HPV80 (3), HPV208 (1), HPV65 (1); 24 jäid alla tüpiseerimisläve.

Sisuliselt: "engineered" ei ole omaette bioloogiline klass, vaid metaandmete lipp, mis eristab
laboris töödeldud rakke naiivsest kliinilisest koest. Enamik neist on rakuliinid või
immortaliseeritud primaarrakud, mistõttu nende kõrge HPV-positiivsus on sama kontaminatsiooni- ja
konstrukti-signaal, mis rakuliinideski.

---

## 2. Kas me teame HPV+ rakuliinide nimesid?

Jah. Vt `tableD_hpv_positive_cell_lines.tsv`. Kokkuvõte (rakuliini tier, HPV+ 93 raamatukogu,
neist 58 tüpiseeritud):

| Rakuliin | Uuring(ud) | Tüüp | Raamatukogusid | Katvus (breadth) | Keskm. sügavus | E1 | E2 | E6 | E7 | L1 |
|:---|:---|:---|---:|:---|:---|---:|---:|---:|---:|---:|
| TIGK | SRP516327, SRP563552 | HPV16 | 18 | 0,10 kuni 0,11 | 45 kuni 514x | 0 | 0 | 59 294 | 34 385 | 0 |
| HK2 | DRP004445 | HPV16 | 6 | 0,10 kuni 0,11 | 52 kuni 242x | 0 | 0 | 9 951 | 14 532 | 0 |
| RWPE-1 | SRP334792 | HPV18 | 6 | 0,15 kuni 0,16 | 109 kuni 183x | 3 304 | 0 | 5 808 | 10 118 | 0 |
| HeLa | SRP418479 | HPV18 | 6 | 0,55 kuni 0,62 | 180 kuni 293x | 122 330 | 109 | 53 670 | 47 573 | 310 |
| SCL-1 | SRP618128, SRP642176 | HPV18 | 12 | 0,56 kuni 0,66 | 41 kuni 238x | 24 646 | 180 | 3 834 | 11 133 | 556 |
| MNT-1 | DRP013144 | HPV18 | 3 | 0,66 | 432 kuni 1319x | 67 752 | 64 | 8 742 | 17 647 | 1 719 |
| HSC-1 | SRP362308 | HPV18 | 2 | 0,35 kuni 0,41 | 3x | 100 | 0 | 3 | 27 | 15 |
| A375 | SRP554912 | HPV18 | 1 | 0,31 | 3x | 51 | 0 | 3 | 16 | 0 |
| A431 | SRP493459 | HPV65 | 1 | 0,74 | 3x | 26 | 5 | 16 | 6 | 7 |

Katvusprofiil eristab kolm selgelt erinevat olukorda, ja see on tõenäoliselt kõige kasulikum osa
viroloogile:

**(a) E6/E7-ainult konstrukt** (breadth 0,10 kuni 0,16, E1 ja E2 ja L1 ja L2 kõik null, sügavus sadu
kordi). 0,105 x 7906 bp on umbes 830 bp, mis vastab täpselt HPV16 E6 + E7 pikkusele.
- **RWPE-1 / HPV18**: teadaolev fakt, ATCC kirjeldab liini kui HPV-18-ga immortaliseeritud
  eesnäärme epiteeli. Meie andmed kinnitavad seda sõltumatult.
- **HK2 / HPV16**: metaandmetes on see kirjas otsesõnu (`cell_type=human papillomavirus 16
  (HPV-16) transformed HK2`).
- **TIGK / HPV16**: see on ootamatu. Kirjanduse järgi on TIGK (telomeraasiga immortaliseeritud
  igeme keratinotsüüdid) tehtud bmi1/hTERT-iga just selleks, et **vältida** HPV E6/E7
  immortaliseerimist. Meie näeme 18 raamatukogus kahest sõltumatust uuringust puhast HPV16 E6+E7
  ekspressiooni (mediaan RPM 255 ja 657). Sisekontroll on tugev: uuringus SRP563552 on samas
  partiis ka HaCaT proovid, ja need on HPV-negatiivsed (mediaan RPM 2,6, ükski ei tüpiseeru).
  See välistab jooksuülese ristsaaste kui ainsa seletuse. Kas kasutatud TIGK alamliin on
  HPV16 E6/E7-transdutseeritud, või on liinivaru saastunud, seda meie andmed ei ütle. **Enne
  publitseerimist vajaks kinnitust** (PCR / STR-tüpiseerimine või uuringute autorite küsitlemine).

**(b) HeLa-taoline osaline genoom** (breadth 0,55 kuni 0,66, E1 domineerib, E2 praktiliselt null,
L1 madal, L2 null). See on HeLa integrandi allkiri: E2 on integratsioonil katkenud, L1/L2 ala suures
osas kadunud.
- HeLa ise annab täpselt selle profiili.
- **SCL-1** (kutaanne SCC liin) ja **MNT-1** (melanoomiliin) annavad sama profiili sama breadth-akna
  sees. Kumbagi ei leidnud ma ICLAC-i valesti identifitseeritud liinide registrist ega kirjandusest
  HPV18-positiivsena. Kõige ökonoomsem seletus on **HeLa ristsaaste**, mis on täpselt see mehhanism,
  mida Cantalupo jt (2015) kirjeldasid TCGA HPV18 valepositiivsete puhul (juba tsiteeritud
  käsikirjas). Melanoomiliinis HPV18 ei ole bioloogiliselt ootuspärane.

**(c) Jälg-tase** (sügavus 3 kuni 4x): HSC-1, A375, A431 (HPV65). Need on liiga madalad, et
konstrukti ja saaste vahel vahet teha. Ma ei paneks neid tabelisse ilma reservatsioonita.

Eraldi: **SRR25439062** on tegelikult nahatüügas (`Skin_wart_18358`, tissue=Skin, HPV2, breadth 0,69,
sügavus 479x, produktiivne), aga see sattus rakuliini tier'i, sest SRA `cell_line` väli sisaldab
väärtust `L11` (ilmselt proovi-, mitte rakuliini tunnus). Vt "Leiud" allpool.

---

## 3. Topeltnakkus positiivsetes kontrollides: millised kombinatsioonid?

Vt `tableC_coinfections.tsv`. Tüpiseeritud produktiivsetest kontrollidest (38 raamatukogu) on
17 (45%) mitme määratud tüübiga. **Aga raamatukogude arv üle hindab sõltumatust:** need 17
raamatukogu taandub 13 unikaalsele lugemisprofiilile ja umbes 7 sõltumatule bioloogilisele proovile.
Põhjus on selgitatud punktis 5: uuringud SRP314958 ja SRP314961 on sama viie patsiendi materjal
(RNA-seq + MeRIP-seq input + IP), ja SRP219878 on 12 tüügast, igaüks kahe jooksuga.

Sõltumatud kombinatsioonid:

| Proov | Kombinatsioon | Perekonnad | Kommentaar |
|:---|:---|:---|:---|
| CA_1 (condyloma) | HPV11 + HPV90 + HPV44 + HPV42 | alpha x4 | kõik neli breadth > 0,87, sügavus 27x kuni 1963x. Kõige kindlam neljikinfektsioon. |
| CA_3 (condyloma) | HPV91 + HPV52 + HPV11 (+HPV42) | alpha x3 (x4) | HPV91 ja HPV52 mõlemad breadth > 0,94, sügavus 83x ja 1114x. HPV11/HPV42 on läve peal. |
| CA_2 (condyloma) | HPV6 + HPV11 | alpha x2 | IP-fraktsioonis muutub HPV6 + HPV42 |
| CA_4, CA_5 (condyloma) | HPV11 + HPV42 | alpha x2 | HPV42 breadth vaid 0,29 kuni 0,59, sügavus 2 kuni 3x |
| 4W (tüügas) | HPV78 + HPV57 | alpha x2 | HPV78 breadth 1,00 sügavus 13 780x; HPV57 breadth 0,66 sügavus 4,6x |
| 60W (tüügas) | HPV27 + HPV-mKN1 | **alpha + gamma** | ainus perekondadeülene kombinatsioon; HPV27 breadth 1,00 sügavus 4998x, HPV-mKN1 breadth 0,48 sügavus 2,2x |

Tähelepanekud, mis on diskussiooni jaoks olulised:

- **Kõik kindlad kombinatsioonid on alpha + alpha.** Ainus perekondadeülene paar (HPV27 + HPV-mKN1,
  alpha + gamma) on nõrga partneriga.
- **Domineeriv tüüp on alati ülekaalukas.** Sekundaarne tüüp on tüüpiliselt 2 kuni 3 suurusjärku
  madalama sügavusega. Ainult CA_1 ja CA_3 puhul on kaks partnerit korraga kõrge katvuse ja
  sügavusega, st ainult neid ma nimetaksin "tõelisteks topeltnakkusteks" ilma reservatsioonita.
- **Läve lähedal olevad kombinatsioonid ei ole stabiilsed.** CA_2 annab RNA-seq ja MeRIP input
  fraktsioonis HPV6 + HPV11, aga IP fraktsioonis HPV6 + HPV42. Sama bioloogiline proov, erinev
  kombinatsioon. See on hoiatus "kokkusobimatute kombinatsioonide" tõlgendamisel: madala katvusega
  partner võib olla ristkaardistus lähedaselt sugulastüübilt (HPV6, HPV11, HPV44, HPV74 on kõik
  alpha-10).
- Rakuliinides ei olnud ühtki ko-infektsiooni (0/58), neoplaasias 3/9.

Nii et jah, saame kombinatsioonid välja tuua, aga ma soovitan need esitada bioloogilise proovi, mitte
raamatukogu tasemel, ja märkida katvus/sügavus iga partneri kohta, et lugeja näeks, kumb on kindel.

---

## 4. Tabel iga üle läve HPV kohta: tüüp / lugemite arv / lokalisatsioon genoomis

Osaliselt on see juba olemas, osaliselt vajab ühte väikest töövoo parandust.

**Olemas kohe:**
- `tableA_assigned_type_calls.tsv`: kõik 223 üle-läve tüübikutset 163 raamatukogus. Veerud: proov,
  tier, uuring, diagnoos, HPV tüüp, perekond, joondunud lugemid, referentsi pikkus, kaetud
  aluspaarid, katvuse laius, keskmine sügavus, RPM.
- `tableB_type_reads_orf.tsv`: 125 raamatukogu, millel on **täpselt üks** määratud tüüp, koos
  ORF-tasemel lugemite arvuga (E1, E2, E5, E6, E7, L1, L2) ja produktiivsuse kutsega. Nendel
  proovidel on ORF-profiil üheselt selle tüübi oma.

**Vajab töövoo parandust (38 raamatukogu, mis on mitme tüübiga):** `featureCounts` jookseb praegu
võtmega `-g gene_name`, mis liidab kõikide referentside samanimelised ORF-id üheks meta-tunnuseks.
Seetõttu on `transcript_classes.tsv` ORF-loendus **kogu proovi kohta, mitte tüübi kohta**. Ko-infektsiooniga
proovis ei saa neid praegu tüüpide vahel jagada.

Parandus on üherealine: GFF-is on juba `ID={referents}_{ORF}` (nt `HPV11REF_E1`), seega piisab
`-g ID` asemel `-g gene_name`, või `-f`, mis annab tunnuse tasemel väljundi koos veergudega
`Chr / Start / End / Strand / Length`. See annaks **korraga** nii ORF-i kui ka nukleotiidide
numeratsiooni, täpselt nagu küsiti. Vt `modules/local/transcript_classify/main.nf` ja
`bin/build_hpv_refs.sh` (GFF-i kirjutamise plokk).

**Nukleotiiditasemel positsioonid** ei ole praegu üldse salvestatud (ainult `covered_bases` ja
`coverage_breadth`). BAM-failid on alles ainult 90 proovil (`results_full_v2/chunks/*/alignments/`);
ülejäänud kustutati koos chunk'i töökataloogidega. Täieliku tabeli saamiseks tuleks tüpiseerimissamm
uuesti läbi lasta ainult HPV+ raamatukogudel (488 tk), mis on odav.

### Seotud viga, mis tuleb igal juhul parandada

`bin/build_hpv_refs.sh` filtreerib PaVE tunnuseid täpse nimega:
`GENES_OF_INTEREST = {"E1","E2","E4","E5","E6","E7","L1","L2"}`.

PaVE aga **ei nimeta neid kunagi nii**:
- E4 on alati `E1^E4` (splaissitud produkt)
- E5 on `E5_ALPHA` / `E5_BETA` / `E5_GAMMA` / `E5_DELTA`

Kontrollisin PaVE API-st (HPV11, HPV16, HPV18, HPV2, HPV12, HPV5, HPV96): mitte üheski ei ole
tunnust nimega `E4` ega `E5`.

Tagajärg andmestikus: E4 rida ei esine **üheski** 488-st `transcript_classes.tsv` failist, ja E5 on
**kõigis 488-s täpselt 0**. Tegelik ORF-profiil on seega ainult E1, E2, E6, E7, L1, L2.

Mõju:
- **Produktiivsuse kutse on puutumata.** See on väravatud ainult L1-le (`late_transcript_min_reads`),
  nii et kõik käsikirja põhijäreldused kehtivad.
- **`EARLY_TOTAL` on alahinnatud**: see on tegelikult E1+E2+E6+E7, mitte E1+E2+E4+E5+E6+E7, nagu
  Methods väidab. E1^E4 on produktiivse infektsiooni kõige rohkearvulisem transkript, nii et selle
  puudumine on sisuliselt märkimisväärne, kui keegi loeb `EARLY_TOTAL`-i kui "kõik varajased".
- Methods'i lauset "aggregated into an early total (E1, E2, E4, E5, E6, E7)" tuleb parandada või,
  parem, viga ära parandada (lisada `E1^E4` ja `E5_*` mustrid) ja E4-ga uuesti jooksutada.

---

## 5. Miks on normaalses nahas kõige sagedasemad HPV11 ja HPV42?

Vastus: **need ei ole normaalne nahk.** See on täielikult ühe uuringu artefakt ja seda ei tohiks
bioloogiliselt tõlgendada.

Kõik 14 HPV11 ja kõik 13 HPV42 kutset "unselected skin" tier'is tulevad kahest SRA uuringust:
**SRP314958** ja **SRP314961**. Sama kehtib HPV74 (5/5), HPV6 (3/3), HPV90 (3/3) ja HPV91 (1/1) kohta.

Mis need proovid on:

- Proovide pealkirjad on `CON_1` kuni `CON_5` ja `CA_1` kuni `CA_5`. `CA` = condyloma acuminatum
  (lesioon), `CON` = sama patsiendi kõrvalkude.
- `CON` proovide metaandmed: `source_name=Skin or mucosa`, `tissue=adjcent mucocal tissues`
  (autorite kirjaviis). See tähendab **kondülooomi kõrval asuv normaalne limaskest**, mitte nahk.
- Meie koekategooria määraja tegi võtmesõnasobitust ja nägi stringis "Skin or mucosa" sõna "Skin",
  seega määras `nahk`. Diagnoos jäi "unspecified", nii et tier'iks sai `unselected_skin`.
- Need proovid tulid **suunatud lesiooni-päringust** (`cohort = lesion`), mitte erapooletust sõelast.
- SRP314958 on RNA-seq; SRP314961 on **MeRIP-seq sama materjali pealt** (input + IP fraktsioonid).
  MeRIP `input` ja RNA-seq annavad identsed lugemiarvud. Seega 5 doonorit esinevad 15 raamatukoguna.

Bioloogiliselt on tulemus täiesti ootuspärane ja isegi huvitav, aga hoopis teise pealkirja all:
**kondülooomi kõrval asuv normaalne limaskest kannab produktiivset HPV11 ja HPV42 infektsiooni.**
Limaskesta tüübid limaskestas. Huvitav detail: lesioonis domineerib HPV11 (sügavus 1963x kuni
12 564x), kõrvalkoes aga HPV42 (breadth kuni 0,999, sügavus 572x kuni 3816x), mis viitab
väljanakkusele (field infection). See on 5 doonoriga ühest uuringust, nii et rohkem kui hüpoteesi
sellest ei saa.

Praktiline järeldus: neid 15 raamatukogu ei tohiks olla "unselected skin" tier'is. Mõju
põhitulemusele on allpool, ja see on oluline.

---

## 6. Joonis alpha/beta/gamma/nu levimusest (Tabel S1 ja S5 juurde)

Tehtud: `figure_genus_prevalence.png` (skript `figure_genus.R`). Kolm rida (määratud tüüp = S1,
jälg-tüüp = S5, ükskõik kumb), viis tier'i, log-skaala, Beta-Binoom posterioorne keskmine 95%
usaldusväärsuse vahemikuga, iga punkti kohal raamatukogude arv. Nullid on märgitud halli ristiga
teljel, mitte "väikese tulbana", sest nulli posterioorne keskmine ei ole null ja tulbana loeks
lugeja seda olemasolevaks signaaliks.

Peamised arvud ("ükskõik kumb tase", % tier'i sees):

| Tier | Alpha | Beta | Gamma | Nu |
|:---|:---|:---|:---|:---|
| Produktiivsed kontrollid (N=43) | **87,5** [76,4 kuni 95,4] | 0 | 10,2 [3,2 kuni 20,6] | 0 |
| Rakuliinid (N=692) | **12,2** [9,9 kuni 14,7] | 0,22 | 0,36 | 0 |
| Neoplaasia (N=442) | 1,92 [0,86 kuni 3,39] | **2,60** [1,33 kuni 4,27] | 0,56 | 0 |
| Valimata nahk (N=7358) | 0,56 | 0,40 | 0,32 | 0,03 |
| Muu (N=2747) | 1,40 | 0,49 | 0,35 | 0 |

Joonise sõnum on väga puhas ja toetab käsikirja: **neoplaasia on ainus tier, kus beta ületab
alpha't**; rakuliinid on praktiliselt puhtalt alpha (HPV16/18); produktiivsed kontrollid on alpha
pluss gamma; nahas on kolm perekonda tasakaalus madalal tasemel. Nu esineb ainult nahas (2
raamatukogu, jälg-tasemel).

Hoiatus: rea "valimata nahk" alpha-osa (0,56%) on üles puhutud punktis 5 ja "Leiud" all kirjeldatud
lekke tõttu. Puhastatud andmetel on nahk sisuliselt beta/gamma-dominantne, mis on **täpselt see
kommensaalne pilt, mida ootaks**, ja see teeb joonise sõnumi tugevamaks, mitte nõrgemaks.

---

## 7. Miks on Tabelis S5 produktiivsete kontrollide veerg 0?

See on konstruktsiooni järgi nii, ja tegelikult on see meetodi positiivne kontroll, mitte puuduv
signaal. Aga tabeli praegune vorming ei ütle seda välja, ja üks arvutus on lisaks eksitav.

**Miks null.** Tabel S5 on *alla läve* jäänud ("jälg") tüübikataloog. Sinna kandideerivad ainult
need HPV+ raamatukogud, mis **ei saanud** määratud tüüpi, st ebaõnnestusid kriteeriumis
breadth >= 0,10 JA sügavus >= 2. Produktiivsed kontrollid läbivad selle kriteeriumi peaaegu alati,
seega nad on Tabelis S1 ja **definitsiooni järgi puuduvad Tabelist S5**.

Arvud:

| Tier | HPV+ | Sai määratud tüübi (S1) | Jäi alla läve (S5-kõlblik) | S5-s lahendus |
|:---|---:|---:|---:|---:|
| Valimata nahk | 234 | 37 | 197 | 56 |
| Rakuliin | 93 | 58 | 35 | 29 |
| Muu | 92 | 21 | 71 | 39 |
| **Produktiivsed kontrollid** | **43** | **38** | **5** | **2** |
| Neoplaasia | 26 | 9 | 17 | 11 |

Nendest 5 alla läve jäänud kontrollist ühel puudub HPV joondus üldse (SRR26678497), kaks kukuvad
läbi jälg-kvaliteedifiltrist (dominantsus < 0,5), ja alles jääb kaks: SRR26678498 (HPV11) ja
SRR26678500 (HPV7). Kuna tabelis kuvatakse ainult tüübid, millel on kokku >= 2 kutset, kaovad need
kaks üksikut enamikust ridadest ära. Seega null.

**Soovitus 1:** kirjuta see allkirja sisse otsesõnu, näiteks: "Productive controls contribute almost
nothing to this table because 38 of 43 productive-control libraries pass the assignment threshold and
therefore appear in Table S1, not here. Their near-absence from the trace tier is a positive control
on the assignment criterion." Ja pane veeru päisesse mõlemad arvud: `Productive controls
(N = 43; 5 trace-eligible)`.

**Soovitus 2 (see on tegelik viga).** Tabeli S5 lahtrites olev "levimus tier'i sees" on arvutatud
**kogu tier'i N kohta**, mitte S5-kõlblike raamatukogude kohta. See on vale nimetaja, sest tier'i
liikmed, kes said S1-s tüübi, ei saanud definitsiooni järgi S5-s esineda. Näiteks HPV16 jälg-kutseid
rakuliinides on 27, ja tabel jagab need 692-ga (3,97%), kuigi jälg-kõlblikke rakuliini raamatukogusid
oli ainult 35 (27/35 = 77%). Erinevus on ligi 20-kordne. Sama kehtib kõigi S5 lahtrite kohta.
Kas tuleks kasutada jälg-kõlblikku nimetajat, või jätta S5-st levimus üldse välja ja näidata ainult
loendeid. Vt `manuscript/_helpers.R`, funktsioon `build_trace_table()` ja `.hdr_N()`.

---

# Leiud, mis nõuavad käsikirjas parandust

Küsimustele vastates tulid välja neli asja. Esimene neist puudutab üht käsikirja väidet otseselt.

## A. "Unselected skin" tier on saastunud kolme uuringurühmaga, ja see kannab M2b tulemust

Praegu on `unselected_skin` tier'is 151 raamatukogu, mis tulid **suunatud lesiooni-päringust**, mitte
erapooletust sõelast. Kolm rühma neist on selgelt valesti paigutatud, ja nad annavad **31 detektsiooni
38-st**, mis ületavad kontrollidega kalibreeritud tuvastuspõranda (13,07 RPM):

1. **SRP314958 + SRP314961**, 15 raamatukogu = 5 doonorit: kondülooomi kõrval asuv limaskest
   (punkt 5). Kõik 15 üle põranda.
2. **SRP394058**, 10 raamatukogu: need on **epidermodysplasia verruciformis'e patsientide
   nahabiopsiad**. `biomaterial_provider=Dr. Vahidnezhad`; uuring on avaldatud kui
   "Whole transcriptome-based skin virome profiling in typical epidermodysplasia verruciformis
   reveals alpha-, beta-, and gamma-HPV infections" (JCI Insight). Neis on produktiivne beta-HPV
   (HPV14 breadth 0,993 sügavus 43 539x, L1 = 334 566 lugemit; ka HPV5, HPV47). Neli neist on
   RPM 34 222 kuni 47 982, st **naha tier'i neli kõrgeimat väärtust on EV biopsiad**. 8 üle põranda.
   Käsikiri väidab praegu (@fig-rpm allkiri): "Epidermodysplasia verruciformis (EV) libraries were
   metadata-flagged as cultures and fall within the cell-line tier." See ei pea paika: kanooniline
   EV biopsiate kohort on `unselected_skin`-is.
3. **SRP560153**, 6 raamatukogu: ENA järgi "Effect of Top2a depletion on the transcriptional profile
   of **HPV31-positive HaCaT (HaCaT31)** cells". Need on HPV31+ rakuliinid shTOP2A knockdown'iga.
   Rakuliini lipp ei tabanud neid, sest metaandmetes on ainult `shTOP2A-1 rep2` ja `skin`, string
   "HaCaT31" ei jõua sõelutud väljadeni. (Kõrvalmärkusena: `cell_line_patterns.py` piirimuster
   `(?![A-Za-z0-9])` ei sobituks stringiga "HaCaT31" niikuinii, sest `hacat` järel on number.)
   Kõik 6 üle põranda.

Tundlikkusanalüüs (tuvastuspõrand 13,07 RPM, sama mudel M2b, log10 sügavus kovariaadina):

| Stsenaarium | Nahk üle põranda | Neoplaasia üle põranda | M2b: neoplaasia vs nahk |
|:---|:---|:---|:---|
| A. nagu käsikirjas | 38/7358 = 0,52% | 6/442 = 1,36% | OR 1,96 [0,81 kuni 4,72], P = 0,14 |
| B. A miinus kondülooomi kõrvalkude | 23/7343 = 0,31% | 6/442 = 1,36% | OR 3,50 [1,38 kuni 8,88], P = 0,0085 |
| C. B miinus EV biopsiad | 15/7333 = 0,20% | 6/442 = 1,36% | OR 6,20 [2,28 kuni 16,9], P = 0,0004 |
| D. C miinus HaCaT31 | 9/7327 = 0,12% | 6/442 = 1,36% | OR 10,0 [3,33 kuni 30,3], P < 0,0001 |
| E. D miinus CIN juhud neoplaasiast | 9/7327 = 0,12% | 3/439 = 0,68% | OR 5,25 [1,33 kuni 20,7], P = 0,018 |

Mida see tähendab:

- **Põhijäreldus ei muutu.** Produktiivsed kontrollid on endiselt 40/43 = 93% üle põranda, st
  ligi kaks suurusjärku üle neoplaasia. M3 produktiivsuse kontrast (L1-põhine) ei muutu üldse, sest
  ükski neist raamatukogudest ei ole kontrollide ega neoplaasia tier'is.
- **Aga üks konkreetne lause käsikirjas läheb ümber.** Praegu (Results, teine alajaotus):
  "neoplasia did not significantly differ from unselected skin after depth adjustment (OR 1,96,
  P = 0,14)". Puhastatud andmetel on neoplaasia **oluliselt** sagedamini üle põranda kui puhas nahk
  (OR 5 kuni 10, P < 0,02 igas variandis). See ei lammuta hit-and-run narratiivi, vaid tugevdab
  seda: neoplaasias on reaalne madala taseme HPV liig puhta naha suhtes, aga see ei jõua kunagi
  produktiivse koormuseni. Bayesi mudel oli seda muide juba vihjanud (P(OR > 1) = 0,9 M2b jaoks).
  Tegelikult muutub lugu sidusamaks: M2a (jälg, "any read") liig ja M2b (põrand) liig osutavad
  samas suunas.
- Praegusel kujul on nahk tehtud "müraseks" lesioonimaterjaliga, mis maskeerib neoplaasia liia.

Soovitus: paranda tier'i määramine (EV biopsiad -> EV/produktiivne kontroll; HaCaT31 -> rakuliin;
kondülooomi kõrvalkude -> eraldi "peri-lesional mucosa" tier või välja) ja jooksuta M2 mudelid uuesti.
Kui EV biopsiad tõsta produktiivsete kontrollide hulka, siis kontrollid on 48/53 = 90,6% üle põranda,
st kontrolli-tier jääb sisuliselt samaks, aga muutub ausamaks (praegu on EV tier tühi).

## B. E4 ja E5 ei ole kunagi loendatud

Vt punkt 4. `EARLY_TOTAL` on tegelikult E1+E2+E6+E7. Produktiivsuse kutse (L1) on puutumata, seega
järeldused kehtivad, aga Methods'i lause on vale ja E1^E4 puudumine on sisuliselt märkimisväärne.

## C. Üks ehtne tüügas on valesti märgitud rakuliiniks

**SRR25439062** (`Skin_wart_18358`, uuring SRP451793): tissue=Skin, diagnoos "control; wart",
HPV2 breadth 0,69, sügavus 479x, E2 = 31 761 lugemit, L1 = 14, produktiivne. See on ilmselgelt
ehtne produktiivne nahatüügas. Ta sattus `cell_line` tier'i, sest SRA `cell_line` väli sisaldab
väärtust `L11`, mis on ilmselt proovi- või liinitähis, mitte tunnustatud rakuliin.
Suund on konservatiivne (kaotame ühe positiivse kontrolli), aga parandada tasuks:
lipp ei tohiks tulla `cell_line` välja pelgast täidetusest, vaid selle sisu peaks mustriga sobituma.

## D. Pseudorepliikatsioon: raamatukogud ei ole sõltumatud

Produktiivsete kontrollide tier'is on 43 raamatukogu, aga umbes **19 sõltumatut bioloogilist proovi**:

- SRP219878: 12 tüügast, igaüks kahe jooksuga (24 raamatukogu, 12 SRX-i)
- SRP314958 + SRP314961: 5 kondülooomi doonorit x 3 raamatukogu (RNA-seq, MeRIP input, MeRIP IP)
- SRP470324: 2 doonorit x 2 tingimust (scRNA-seq)

Sama kehtib naha tier'is CON proovide kohta (5 doonorit, 15 raamatukogu). Praegused Wilsoni ja
Beta-Binoomi usaldusvahemikud eeldavad raamatukogude sõltumatust, seega on nad **liiga kitsad**.
Miinimum: mainida seda Methods'is / limitatsioonides. Parem: koondada SRX või doonori tasemele, või
kasutada uuringu juhuslikku efekti (glmer / rstanarm `(1 | study)`).

---

---

# Tier'i parandus ja M2 mudelite uuesti jooksutamine

Punktis A kirjeldatud leke on nüüd parandatud töövoo tasemel. Vana väljund
(`results_reframe_csc/`) on puutumata, uus on `results_reframe_fixed/`. Käsikirja teksti ma ei
muutnud, sest see nõuab sisulist otsust (vt allpool).

## Mis muutus koodis

- **`assets/tier_overrides.csv`** (uus). Kureeritud paranduste tabel, iga rida koos tõendiga.
  Veerud: `scope` (sample või study), `id`, `field` (`tier` või `is_cell_line`), `value`, `reason`.
  18 reeglit:
  - `SRP560153` -> `is_cell_line = TRUE` (HPV31+ HaCaT31 rakud; ENA uuringu pealkiri)
  - `SRR25439062` -> `is_cell_line = FALSE` (ehtne nahatüügas; SRA `cell_line` väli sisaldab
    proovitähist "L11")
  - `SRP394058` -> `tier = EV` (EV patsientide nahabiopsiad; SRA-s puudub diagnoosi võtmesõna)
  - 15 x `CON_*` proovi (`SRP314958`, `SRP314961`) -> `tier = peri_lesional`
- **`bin/hpv_rpm_contrast.R`**. Uus võti `--tier-overrides` (vaikimisi
  `assets/tier_overrides.csv`). `is_cell_line` parandused rakenduvad **enne** tier'i reegleid
  (nad toidavad ülimuslikkuse ahelat), `tier` parandused **pärast** (nad on viimane sõna).
  Proovitasand võidab uuringutasandi. Lisatud uus tier `peri_lesional`. `stopifnot` valvab, et
  ükski override ei too sisse tundmatut tier'i.
- **`manuscript/collab_answers/rerun_m2.R`** (uus). Sobitab M2 mõlemale tabelile ja trükib
  kõrvuti võrdluse. Väljund: `m2_rerun_output.txt`.

Kontroll: `bin/hpv_rpm_contrast.R --tier-overrides ""` taastoodab `results_reframe_csc`
täpselt (11 282 rida, identsed tier'id, RPM erinevus 0). Parandus ei muuda midagi peale nende
31 raamatukogu.

## Tier'ide koosseis enne ja pärast

| Tier | Enne (n) | Pärast (n) | Üle põranda enne | Üle põranda pärast |
|:---|---:|---:|---:|---:|
| Valimata nahk | 7358 | **7327** | 0,5% | **0,1%** |
| Rakuliin | 692 | 697 | 10,3% | 10,9% |
| Produktiivsed kontrollid | 43 | 44 | 93,0% | 93,2% |
| EV | 0 | **10** | n/a | **80,0%** |
| Peri-lesionaalne limaskest | 0 | **15** | n/a | **100%** |
| Neoplaasia | 442 | 442 | 1,4% | 1,4% |
| Muu | 2747 | 2747 | 1,0% | 1,0% |

Tuvastuspõrand nihkus 13,07 -> 13,22 RPM (produktiivsete kontrollide hulka lisandus üks ehtne
tüügas). Uus joonis 1 on `manuscript/figure1.png` (renderdatud parandatud tier'itega).

## M2 mudelite tulemused

| Kontrast | Enne (avaldatud) | Pärast (parandatud) |
|:---|:---|:---|
| M2a neoplaasia vs nahk, ükskõik milline HPV lugem | OR 1,66 [1,36 kuni 2,03], P = 7,9e-07 | OR 1,69 [1,38 kuni 2,07], P = 3,3e-07 |
| M2a sekveneerimissügavus (10x lugemite kohta) | OR 4,65 [4,08 kuni 5,30] | OR 4,59 [4,03 kuni 5,23] |
| **M2b neoplaasia vs nahk, tuvastuspõrandal** | **OR 1,96 [0,81 kuni 4,72], P = 0,14** | **OR 10,04 [3,33 kuni 30,25], P = 4,1e-05** |
| M2b sama, CIN juhud neoplaasiast eemaldatud | OR 0,98 [0,30 kuni 3,23], P = 0,98 | OR 5,25 [1,33 kuni 20,71], P = 0,018 |
| M2c kontrollid vs nahk, põrandal (Fisher) | OR 2515 [754 kuni 16383], P = 4,9e-81 | OR 8295 [2409 kuni väga suur], P = 3,3e-96 |
| M2d sügavus neoplaasia sees, põrandal | OR 1,71 [0,13 kuni 22,58], P = 0,68 | muutumatu |
| EV vs nahk, põrandal (Fisher, eraldatus) | n/a | OR 3093 [508 kuni väga suur], P = 5,3e-21 |
| Peri-lesionaalne vs nahk (Fisher, täielik eraldatus) | n/a | OR lõpmatu [2409 kuni lõpmatu], P = 1,8e-40 |

Mudelid 1 (positiivsus proovi klassi kaupa) ja 3 (L1-põhine produktiivsus) ei loe tier'i tabelit,
seega **need ei muutu üldse**.

## Mida see tähendab

**Põhijäreldus püsib.** Produktiivsed kontrollid on 93% üle põranda, neoplaasia 1,4%, st ligi kaks
suurusjärku vahet. M2c kontrast tugevneb (OR 2515 -> 8295). M3 produktiivsuse kontrast on puutumata.
"Trace, mitte produktiivne" lugu jääb kehtima.

**Aga kaks käsikirja väidet lähevad ümber.**

1. *"Neoplasia did not significantly differ from unselected skin after depth adjustment."*
   Nüüd erineb: OR 10,0, P = 4,1e-05. Ka pärast CIN juhtude eemaldamist OR 5,25, P = 0,018.
   Vana tulemus tekkis sellest, et nahk oli lesioonimaterjaliga "mürastatud": 31 raamatukogu
   38-st üle põranda olevast naha raamatukogust olid tegelikult EV biopsiad, kondülooomi
   kõrvalkude või HPV31+ rakuliin.

   See **ei lammuta** hit-and-run mudelit, vaid muudab loo sidusamaks. Enne oli veidi ebamugav, et
   M2a (jälg-tase, OR 1,66, P < 1e-6) näitas neoplaasias liia, aga M2b (põrand) mitte. Nüüd
   osutavad mõlemad samas suunas: neoplaasias on reaalne madala taseme HPV liig puhta naha suhtes
   (OR 1,7 jälg-tasemel, OR 5 kuni 10 põranda tasemel), ja see jääb ikkagi kaks suurusjärku alla
   produktiivsete kontrollide taseme. Bayesi mudel oli seda muide juba vihjanud: M2b posterioorne
   P(OR > 1) oli 0,9.

   Uus sõnastus võiks olla umbes: neoplaasia kannab HPV-d oluliselt sagedamini kui puhas nahk, aga
   tuvastatav koormus jääb kaugele alla produktiivse infektsiooni taseme, mida näeme tüügastes ja
   EV-s. See on täpselt "residuaalne jalajälg, mitte aktiivne infektsioon".

2. *"Epidermodysplasia verruciformis (EV) libraries were metadata-flagged as cultures and fall
   within the cell-line tier."* (@fig-rpm allkiri.) See ei pea paika. EV on nüüd oma tier,
   n = 10, 80% üle põranda, mediaan RPM 230, ja see on **kõige sügavamalt sekveneeritud tier**
   (mediaan 62,6M lugemit). EV biopsiad kuuluvad joonisele omaette veeruna, kus nad on selge
   positiivne kontroll produktiivsete kontrollide kõrval.

Puudutatud kohad käsikirjas (kõik viitavad `stats$m2b_or` / `stats$m2b_p`):
`_abstract.qmd` (1 koht), `_results.qmd` (2 kohta: peatükk "Skin neoplasia carries only trace HPV"
ja Bayesi peatükk), `_discussion.qmd` (1 koht), pluss @fig-rpm allkirja tier'ide N-id ja EV lause,
ning `stats$uns_pctfloor` (0,5 -> 0,1) kahes kohas.

---

# Käsikirja parandus ja Bayesi mudelite ümbersobitus

Kogu analüüsiahel on nüüd parandatud tier'itel ja käsikiri renderdub puhtalt (HTML + DOCX,
ilma hoiatusteta).

## Mis koodis muutus

- **`manuscript/_paths.R`** (uus). Üks tõe allikas: `REFRAME_DIR <- "results_reframe_fixed"`,
  pluss `apply_cellline_overrides()`, mis viib samad kureeritud `is_cell_line` parandused ka
  mudelisse 1. Varem oleks HaCaT31 olnud tier'i tabelis rakuliin, aga mudelis 1 "kliiniline kude".
- **`results_models.R`**: loeb parandatud tabelit, rakendab rakuliini override'id mudelile 1,
  lisab EV ja peri-lesionaalse tier'i statistikasse, arvutab EV ja peri-lesionaalse Fisheri
  kontrastid (mõlemad on nahaga eraldatud, seega glm ei kõlba), ning ei trüki enam eraldatuse
  korral 16-kohalist usaldusvahemiku ülempiiri (nüüd `∞`).
- **`bayesian_models.R`**: loeb parandatud tabelit, samad override'id, sama mudelistruktuur.
- **`_helpers.R`**: tabelitesse S1 ja S5 lisandusid EV ja peri-lesionaalse veerud; `fo()` oskab
  vormindada lõpmatust; uus `trace_elig_str()` trükib S5 allkirja jaoks jälg-kõlblike arvu tier'i
  kohta.
- **`references.bib`**: lisatud Saeidian jt 2023 (JCI Insight), EV biopsiate kohordi allikas.
- `manuscript/figure1.png` ja `figureS1.png` on uuesti genereeritud parandatud tier'itega.

## Bayesi posterioorid (nõrgalt informatiivsed priorid, Normal(0, 2.5))

| Kontrast | Enne | Pärast |
|:---|:---|:---|
| M1 rakuliin vs kliiniline | 4,4 [3,0 kuni 6,3] | 5,2 [3,6 kuni 7,3] |
| M1 engineered vs kliiniline | 5,5 [3,3 kuni 8,6] | 4,6 [2,8 kuni 7,4] |
| M2a neoplaasia vs nahk, jälg | 1,7 [1,4 kuni 2,1] | 1,7 [1,4 kuni 2,1] |
| **M2b neoplaasia vs nahk, põrand** | **1,9 [0,7 kuni 4,2], P(OR>1) = 0,90** | **8,7 [2,9 kuni 26], P(OR>1) > 0,999** |
| M2b sügavus, põrand | 2,8 [1,4 kuni 5,9], P(OR>1) = 0,9995 | 1,5 [0,5 kuni 5,1], P(OR>1) = 0,74 |
| M2c kontrollid vs nahk, põrand | 1 837 [715 kuni 6 087] | 6 709 [2 330 kuni 23 011] |
| M2c sama, Student-t(3) prior | 2 527 [873 kuni 10 576] | 10 715 [3 250 kuni 46 581] |
| M2d sügavus neoplaasia sees | 1,8 [0,2 kuni 21], P(OR>1) = 0,70 | muutumatu |
| M3 kontrollid vs neoplaasia | 3,7 [1,3 kuni 11], P(OR>1) = 0,991 | muutumatu |

Kõik ahelad koondusid (R-hat ≤ 1,002, divergentseid siirdeid ei olnud). Kaks asja väärivad
tähelepanu. Esiteks, M2b posterior ei sisalda enam nulli: väike, aga tõeline põranda-tasemel liig
neoplaasias on kinnitatud ka Bayesi raamistikus. Teiseks, sügavuse kovariaat **kadus** põranda
tasemel (2,8 -> 1,5, P(OR>1) langes 0,9995-lt 0,74-le). Varem tundus, et sügavam sekveneerimine
tõstab põrandani jõudmise šansse, aga see oli tegelikult lekke artefakt: EV biopsiad ja
peri-lesionaalne kude olid nii kõige sügavamalt sekveneeritud kui ka kõik üle põranda. Kui nad
nahast välja tõsta, kaob sügavuse efekt põranda tasemel täielikult, mis **tugevdab** käsikirja
keskset argumenti, et vahe ei tulene sügavusest.

## Mis tekstis muutus

- **Abstrakt**: "did not differ from unselected skin at productive level" asendatud. Uus sõnastus
  ütleb, et neoplaasia jõuab põrandani 1,4% raamatukogudest vs kontrollide 93,2%, mediaankoormus on
  kuus suurusjärku madalam, ja et sügavusega korrigeerituna on neoplaasia küll puhtast nahast
  kõrgemal, aga ainult jäljena. "four- to five-fold" -> "roughly five-fold", "P < 10⁻¹¹" -> "P < 10⁻⁸".
- **Tulemused, M2 lõik**: ümber kirjutatud kolmeks lõiguks. Uus rõhk on absoluutskaalal (41/44 vs
  6/442 vs 9/7327), sest OR-id on suured ainult seetõttu, et naha baasjoon on nulli lähedal.
  Lisatud, et sügavus ei tõsta põrandani jõudmise šansse ei kovariaadina (OR 1,4, P = 0,59) ega
  neoplaasia sees (OR 1,7, P = 0,68), ja et kõige sügavamad cSCC/AK raamatukogud (327,3M lugemit)
  jäävad allapoole põrandat.
- **Tulemused, CIN lõik**: nüüd öeldakse otse, et 3 kuuest põrandaületajast on limaskesta HPV16
  juhud, ja et nende eemaldamine poolitab efekti (OR 10,0 -> 5,3). Seega on peatekstis raporteeritud
  liig ülemine hinnang.
- **Tulemused, Bayesi lõik**: varem "the only effect the posterior nuances...". Nüüd: posterior
  kahandab efekti, aga jätab selgelt positiivseks.
- **Tulemused, tüübimaastik**: uus tekst kirjeldab, et EV tier kannab täpselt neid beeta-tüüpe, mis
  sündroomi defineerivad (HPV9, 14, 20, 25, 47), peri-lesionaalne limaskest kannab kõrvalasuva
  kondülooomi alfa-repertuaari, ja valimata nahk ei anna enam ühtki tüüpi üle kuvamisläve peale kahe
  HPV18 kontaminatsioonikutse. Naha beeta/gamma reservuaar on tervikuna jälg-tasemel (tabel S5).
- **Joonise 1 allkiri**: uued tier'id ja nende N-id; vale lause EV kohta asendatud korrektsega
  (EV *kultuurid* jäävad rakuliini tier'i, EV tier hoiab 10 patsiendi nahabiopsiat).
- **Joonise S1 allkiri**: mainib nüüd EV ja peri-lesionaalset ning selgitab, et tier joonistatakse
  ainult neis sügavusvahemikes, kus tal on vähemalt kolm raamatukogu.
- **Meetodid**: uus lõik, mis kirjeldab kõiki nelja kureeritud parandust koos tõenditega, ütleb, et
  need liigutavad 31 raamatukogu 11 282-st, et samad `is_cell_line` parandused rakenduvad ka
  mudelile 1, ja et override'ide väljalülitamine taastoodab vana tier'i määramise täpselt.
- **Tabel S1 allkiri**: võrdlus "14 HPV11 kontrollides vs 14 nahas" oli parandatud tiering'us vale;
  nüüd võrreldakse kontrolle peri-lesionaalse limaskestaga (32% vs 91%).
- **Tabel S4 allkiri**: "indistinguishable from unselected skin" eemaldatud.
- **Tabel S5 allkiri**: vastab nüüd otse kaasautori küsimusele. Öeldakse, et S5 on tabeli S1
  **täiend**, mitte laiendus, et iga raamatukogu esineb täpselt ühes neist, ja et 44-st
  produktiivsest kontrollist 39 saab tüübi tabelis S1, jättes S5-le 5, millest 2 lahenduvad. Null on
  seega meetodi positiivne kontroll. Lisatud jälg-kõlblike arv tier'i kohta, et lugeja saaks
  protsente ise ümber arvutada.
- **Diskussioon**: M2b lause ümber kirjutatud; "most deeply sequenced tier" asendatud (EV biopsiad
  on tegelikult sügavaim tier, mediaan 62,6M); "four- to five-fold" -> "roughly five-fold".

## Mida ma endiselt EI teinud

- **Ei muutnud küsimuste 1 kuni 7 tabeleid ega perekonna-joonist** (`tableA` kuni `tableD`,
  `figure_genus_prevalence.png`). Need on tahtlikult **avaldatud** tier'ite peal, sest just nemad
  dokumenteerivad probleemi: näiteks HPV11 ja HPV42 ilmuvad neis veel "unselected skin" all, mis on
  kogu 5. küsimuse vastuse tõendusmaterjal. Parandatud tier'itel jooksutamiseks vaheta
  `make_tables.R` ja `figure_genus.R` sees `results_reframe_csc` -> `results_reframe_fixed`.
- **Ei parandanud E4/E5 viga** (leid B). See nõuab referentside uuesti ehitamist ja HPV+ proovide
  tüpiseerimissammu kordamist, ega mõjuta ühtki praegust järeldust (produktiivsuse kutse on
  L1-väravatud).
- **Ei tegelenud pseudorepliikatsiooniga** (leid D). Usaldusvahemikud eeldavad endiselt
  raamatukogude sõltumatust.

---

## Failid

| Fail | Sisu |
|:---|:---|
| `tableA_assigned_type_calls.tsv` | 223 üle-läve tüübikutset, 163 raamatukogu (küsimus 4) |
| `tableB_type_reads_orf.tsv` | 125 ühe tüübiga raamatukogu, ORF-tasemel lugemid (küsimus 4) |
| `tableC_coinfections.tsv` | 38 ko-infektsiooniga raamatukogu, kombinatsioonid + katvus (küsimus 3) |
| `tableD_hpv_positive_cell_lines.tsv` | HPV+ rakuliinid nimeliselt, tüüp + katvusprofiil (küsimus 2) |
| `figure_genus_prevalence.png` | perekonna levimus tier'i kaupa, S1 ja S5 (küsimus 6) |
| `m2_rerun_output.txt` | M2 mudelite täisväljund, enne ja pärast |
| `make_tables.R`, `figure_genus.R`, `rerun_m2.R` | reprodutseerimiseks, jooksutada repo juurkataloogist |

Väljaspool seda kausta: `assets/tier_overrides.csv` (uus), `bin/hpv_rpm_contrast.R` (muudetud),
`results_reframe_fixed/` (uus väljund).
