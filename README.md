# hse21_H3K4me3_G4_human

## Майнор "Биоинформатика", 2021

#### Выполнил

Серебренников Александр Дмитриевич

Группа БПИ181, МБ-2

Сервер: 92.242.58.92, Порт: 5222

Исходные файлы: human (hg19), DNA structure - 
[G4_seq_Li_K](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3003539), 
гистоновая метка - [H3K4me3](https://www.encodeproject.org/chip-seq-matrix/?type=Experiment&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens&assay_title=Histone%20ChIP-seq&assay_title=Mint-ChIP-seq&status=released),
тип клеток - Н1, BED-файлы с данными экспериментов - 
[ENCFF178RTX](https://www.encodeproject.org/experiments/ENCSR019SQX/) 
и [ENCFF236DPM](https://www.encodeproject.org/experiments/ENCSR003SSR/).

#### Анализ пиков гистоновой метки

Загружаем данные экспериментов в формате BED и оставляем в файлах первые 5 столбцов. 
Геномные координаты в данных файлах уже соответствуют референсу человеческого генома 
hg19, поэтому переводить их в другие не требуется. 

```bash
wget https://www.encodeproject.org/files/ENCFF178RTX/@@download/ENCFF178RTX.bed.gz
wget https://www.encodeproject.org/files/ENCFF236DPM/@@download/ENCFF236DPM.bed.gz
zcat ENCFF178RTX.bed.gz  |  cut -f1-5 > H3K4me3_H1.ENCFF178RTX.hg19.bed
zcat ENCFF236DPM.bed.gz  |  cut -f1-5 > H3K4me3_H1.ENCFF236DPM.hg19.bed
```

Для обработки и анализа геномных данных используется 
[проект](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/tree/main/src) 
на R. Также, здесь и далее передача данных с локального компьютера на сервер и наоборот
осуществляется через коммиты и стягивание изменений из этого репозитория.

[Данный скрипт](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/src/lib.R) - технический, 
нужен для загрузки и подключения необходимых библиотек и установки путей. Он используется другими скриптами.

Далее, с помощью [скрипта](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/src/len_hist.R)
построим гистограмму длин участков для каждого эксперимента. 
На графиках также отражено кол-во пиков в каждом файле.

![len_hist.H3K4me3_H1.ENCFF178RTX.hg19](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/images/len_hist.H3K4me3_H1.ENCFF178RTX.hg19.png)

![len_hist.H3K4me3_H1.ENCFF236DPM.hg19](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/images/len_hist.H3K4me3_H1.ENCFF236DPM.hg19.png)

Теперь, из двух данных файлов с ChIP-seq пиками выкидываем слишком длинные пики (outliers). 
С помощью [этого скрипта](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/src/filter_peaks.R)
оставим в файле H3K4me3_H1.ENCFF178RTX.hg19 все пики короче 4500, 
а в файле H3K4me3_H1.ENCFF236DPM.hg19 оставим все пики короче 5500. Данные значения
были определены по гистограммам и отсортированному содержимому файлов. 
Распределение длин пиков в обновленных файлах:

![filter_peaks.H3K4me3_H1.ENCFF178RTX.hg19.filtered.hist](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/images/filter_peaks.H3K4me3_H1.ENCFF178RTX.hg19.filtered.hist.png)

![filter_peaks.H3K4me3_H1.ENCFF236DPM.hg19.filtered.hist](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/images/filter_peaks.H3K4me3_H1.ENCFF236DPM.hg19.filtered.hist.png)

Затем, с помощью [следующего скрипта](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/src/chip_seeker.R)
и R-библиотеки ChIPseeker проанализируем, где располагаются пики 
гистоновой метки относительно аннотированных генов. 
Получившиеся пай-чарт графики:

![chip_seeker.H3K4me3_H1.ENCFF178RTX.hg19.filtered.plotAnnoPie](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/images/chip_seeker.H3K4me3_H1.ENCFF178RTX.hg19.filtered.plotAnnoPie.png)

![chip_seeker.H3K4me3_H1.ENCFF236DPM.hg19.filtered.plotAnnoPie](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/images/chip_seeker.H3K4me3_H1.ENCFF236DPM.hg19.filtered.plotAnnoPie.png)

Далее на сервере два отфильрованных файла объединяются с помощью утилиты bedtools, команды merge 
(bedtools merge не принимает не отсортированный файл, поэтому мы предварительно его сортируем).

```bash
adserebrennikov@laboratory01:~/ngs/project-git/hse21_H3K4me3_G4_human/data$ cat  *.filtered.bed  |   sort -k1,1 -k2,2n   |   bedtools merge   >  H3K4me3_H1.merge.hg19.bed 
```

Визуализируем исходные два набора ChIP-seq пиков, а также их объединение в геномном браузере.
Ссылка на соответствующую [сессию](https://genome.ucsc.edu/s/SerebrennikovAlexandr/minor_bioinformatics). Ниже приложен скриншот.

![Скриншот_сессии_в_геномном_браузере](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/screenshots/session1.png)

Видно, что bedtools merge сработал корректно, объединенный файл - действительно объединение пиков из двух файлов.

#### Анализ участков вторичной структуры ДНК

Теперь необходимо получить и проанализировать данные вторичной структуры ДНК.

Скачиваем два BED-файла со вторичной структурой ДНК, оставляем первые 5 столбцов и объединяем оба файла аналогично файлам с гистоновыми метками.

```bash
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3003nnn/GSM3003539/suppl/GSM3003539_Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3003nnn/GSM3003539/suppl/GSM3003539_Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed.gz
zcat GSM3003539_Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed.gz | cut -f1-5 > GSM3003539_plus.bed
zcat GSM3003539_Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed.gz | cut -f1-5 > GSM3003539_minus.bed
cat GSM3003539_*.bed | sort -k1,1 -k2,2n | bedtools merge > GSM3003539.merged.bed 
```

Также, аналогично файлам с гистоновыми метками, 
строим распределение длин участков вторичной структуры вместе с
подсчетом количества пиков. Смотрим, где располагаются участки стр-ры ДНК относительно аннотированных генов.

![len_hist.GSM3003539.merged](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/images/len_hist.GSM3003539.merged.png)

![chip_seeker.GSM3003539.merged.plotAnnoPie](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/images/chip_seeker.GSM3003539.merged.plotAnnoPie.png)

#### Анализ пересечений гистоновой метки и структуры ДНК

Далее найдем пересечения между гистоновой меткой и структурами ДНК с помощью bedtools intersect,
найдем количество и распределение длин этих пересечений:

```bash
  bedtools intersect -a GSM3003539.merged.bed -b H3K4me3_H1.merge.hg19.bed > intersect.bed
```

![len_hist.intersect](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/images/len_hist.intersect.png)

Затем визуализируем в геномном браузере исходные участки 
структуры ДНК, а также их пересечения с гистоновой меткой.
Ссылка на соответствующую [сессию](https://genome.ucsc.edu/s/SerebrennikovAlexandr/minor_bioinformatics_2).

Для дальнейшего анализа нам интересны места, где есть 
пересечение между гистоновой меткой и стр-рой ДНК 
(желательно рядом с аннотированным геном). 
Ниже в качестве примера скриншоты и геномные координаты этих мест.

###### chr1:74663447-74663571

![Скриншот_первого_участка_в_геномном_браузере](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/screenshots/session2.png)

###### chr1:75139011-75139100

![Скриншот_второго_участка_в_геномном_браузере](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/screenshots/session3.png)

Ассоциируем полученные пересечения с ближайшими генами c 
помощью [скрипта](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/src/ChIPpeakAnno.R).
В результате его работы получаем 2 файла. 
Первый - [файл ассоциаций пиков с генами](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.intersect_with_G4_seq_Li_K.genes.txt) (14906 строк), 
а также [список уникальных генов](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/data/H3K4me3_H1.intersect_with_G4_seq_Li_K.genes_uniq.txt) (8961 строка). 

GO-анализ для полученных уникальных генов.

![GO_результат](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/screenshots/res1.png)

![GO_результат_таблица](https://github.com/SerebrennikovAlexandr/hse21_H3K4me3_G4_human/blob/main/screenshots/res2.png)
