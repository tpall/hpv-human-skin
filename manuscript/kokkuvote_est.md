---
title: "HPV inimese nahas ja nahakasvajates — tulemuste kokkuvõte"
subtitle: "Avalike RNA-seq andmete taasanalüüs · sisend kaasautoritelt"
date: "5. juuli 2026"
lang: et
---

## Taust ja eesmärk

Enamik inimese papilloomiviirustest (HPV) on naha beeta- ja gammatüüpi viirused, mis püsivad nahal peamiselt kommensaalse, subkliinilise reservuaarina. Nende roll naha lamerakk-kartsinoomis (cSCC) on vaieldav — "hit-and-run" onkogeneesimudeli ja kaitsva immuunsuse mudeli vahel. Analüüsisime süstemaatiliselt avalikke inimese RNA-seq andmeid (nahk ja seotud koed), et kirjeldada, millised HPV tüübid nahas esinevad, kui sageli ja kas nad on transkriptsiooniliselt aktiivsed (produktiivsed).

## Meetod

Ühtne töövoog skriinis iga teegi papilloomiviiruse lugemite suhtes (Kraken2, inimese genoom + HPV referentspaneel), tüpiseeris positiivsed kogu paneeli vastu (STAR + katvuse/sügavuse künnised) ja klassifitseeris infektsiooni produktiivseks hiliste (L1, peamise kapsiidivalgu) transkriptide alusel. Rakuliinid ja muundatud materjal märgistati eraldi; iga teek skooriti sekveneerimissügavusest sõltumatult (lugemid miljoni kohta, RPM). Kokku skriiniti **11 556 teeki** (HPV-positiivseid 2,54%).

## Tulemus 1 — signaal on koondunud mitte-primaarsesse materjali (kontaminatsioon)

HPV-positiivsus oli kliinilises primaarkoes vaid **2,15%**, kuid muundatud teekides **11,1%** (šansisuhe 5,7) ja rakuliinides **9,0%** (šansisuhe 4,5) — 4–5 korda kõrgem (p = 1,4×10⁻¹⁸). Need positiivsed olid valdavalt integreeritud **HPV16/HPV18** — HeLa-tüüpi rakuliinide ja nende hästi dokumenteeritud arhiivikontaminatsiooni allkiri. Avalikes transkriptoomides nähtav "HPV rohkus" on seega suuresti laboriartefakt, mis tuleb enne igasugust levimushinnangut eraldada.

## Tulemus 2 — nahakasvajates on HPV vaid jäljetasemel, mitte produktiivne (põhitulemus)

cSCC/aktiinilise keratoosi (neoplaasia, n = 442) teekides leidus HPV-lugemeid pooltes proovides (~50%), kuid produktiivsel tasemel vaid **1,4%-l** (mediaan ~0,005 RPM) — hoolimata sellest, et see oli **kõige sügavamalt sekveneeritud** kude. Sügavuse järgi sobitatud võrdlus on otsustav: erinevus produktiivsete kontrollide (tüükad/EV: 100% avastatud, mediaan ~5400 RPM) ja naha vahel püsib (Fisheri šansisuhe **2515**, p = 5×10⁻⁸¹), samas kui neoplaasia **ei erine tervest nahast** produktiivsel tasemel (sügavuse suhtes korrigeeritult mitteoluline, p = 0,14). Beeta-HPV cSCC-s on seega **jäljetranskriptsioon**, kooskõlas hit-and-run / UV-kofaktori mudeliga, mitte tüügaste produktiivse koormusega. See laiendab varasemat leidu (Arron 2011) umbes kaks suurusjärku suuremal skaalal ja sügavuse suhtes kontrollituna.

## Tulemus 3 — produktiivne infektsioon on tüügaste ja EV omadus

Produktiivseid infektsioone (L1 ≥ 3 lugemit) leidus rikkalikult tüükastes ja epidermodysplasia verruciformis'es (**64%**) vs neoplaasias (**32%**; šansisuhe 3,8, p = 0,013). Siiski esines üksikuid ehtsaid, mitte-rakuliinist produktiivseid koeproove: **HPV12** cSCC-s (katvus 99,8%, L1 = 1803), **HPV96 + HPV182** koinfektsioon (aktiiniline keratoos/kartsinoom, katvus 98–99%) ja **HPV16**. Populatsioonitasandi "jälje"-järeldus kirjeldab cSCC/AK keskmist, mitte seda, et produktiivset naha-HPV haigust üldse ei esineks.

## Robustsus

Kõik kolm järeldust püsivad Bayesi logistilise regressiooni all nõrgalt informatiivsete prioorite (Normal(0; 2,5)) korral.

## Järeldus

Naha-HPV on **produktiivne tüükastes ja EV-s, kuid transkriptsiooniliselt vaikne tervel nahal ja kartsinoomis**, kus ta esineb vaid jälje, mitteproduktiivse "jalajäljena". Seda fenomenoloogiat peavad seletama nii hit-and-run kui ka kaitsva immuunsuse mudel. Laiemalt pakub kontaminatsiooniteadlik, sügavus-normaliseeritud raamistik skaleeritava viisi avalike transkriptoomide "äravisatud" (mitte-inimese) fraktsiooni kaevandamiseks.
