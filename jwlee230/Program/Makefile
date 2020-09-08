# Variables
DATE := $(shell date "+%y%m%d")
RANDOM := $(shell bash -c 'echo $$RANDOM')
IMAGE_NAME = premature:latest
CPUS = 20
MEMS = 100G
PWD := $(shell pwd)
TOOLS = 
DOCKER = docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) $(IMAGE_NAME)

# Options
VOLUME_OPTS = --volume $(abspath Output):/Output --volume $(abspath Data):/Data --volume $(abspath Python):/Python
RUN_OPTS = --tty --cpus="$(CPUS)" --memory="$(MEMS)"

# General
all:
.PHONY += all

log Output Docker/Tools:
	mkdir $@

# Tools
Docker/Tools/Mothur.zip: | Docker/Tools
	wget "https://github.com/mothur/mothur/releases/download/v1.44.2/Mothur.Ubuntu_18.zip" -O $@
TOOLS += Docker/Tools/Mothur.zip

# Docker
TOOLS += $(wildcard Docker/*)
build.log: Docker/Dockerfile $(TOOLS) | log Output
	rm -fv $@
	docker images | grep $(IMAGE_NAME) && docker rmi $(IMAGE_NAME) || true
	docker build --rm --tag $(IMAGE_NAME) $(<D) | tee $@

build: build.log
.PHONY += build

interactive: build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) --interactive $(IMAGE_NAME) /bin/bash || true
.PHONY += interactive

delete: build.log
	docker rmi $(IMAGE_NAME)
	rm -fv build.log
.PHONY += delete

stop:
	docker rm $(CONTAINER_NAME)

# SGE
tmp.sh: build.log
	echo "make -C $(PWD) latest" > tmp.sh

run: tmp.sh | log Output
	qsub -cwd -l h_vmem=$(MEMS) -m abe -M "230@fumire.moe" -N Premature_$(DATE) -pe smp $(CPUS) -o $(abspath log) -e $(abspath log) $<
.PHONY += run

# Actual
latest: step09
.PHONY += latest

# Step 01 (Make manifest files)
Output/Manifest:
	mkdir -p $@

Output/Manifest/whole.tsv: Python/step01.py $(wildcard Data/RawData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

step01: Output/Manifest/whole.tsv
.PHONY += step01

# Step 02 (De-multiplexing)
Output/Demux:
	mkdir -p $@

Output/Demux/%.qza: Output/Manifest/%.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/%.qzv: Output/Demux/%.qza | build.log
	$(DOCKER) qiime demux summarize --i-data $(addprefix /,$<) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step02: Output/Demux/whole.qza Output/Demux/whole.qzv
.PHONY += step02

# Step 03 (Quality filter)
Output/Filter:
	mkdir -p $@

Output/Filter/%.seq.qza Output/Filter/%.stat.qza: Output/Demux/%.qza | Output/Filter build.log
	$(DOCKER) qiime quality-filter q-score --i-demux $(addprefix /,$<) --o-filtered-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-filter-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step03: Output/Filter/whole.seq.qza
.PHONY += step03

# Step 04 (DADA2 Denoising)
Output/DADA2:
	mkdir -p $@

Output/DADA2/%.table.qza Output/DADA2/%.seq.qza Output/DADA2/%.stat.qza: Output/Demux/%.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f 265 --p-trunc-len-r 265 --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step04: Output/DADA2/whole.seq.qza
.PHONY += step04

# Step 05 (Deblur Denoising)
Output/Deblur:
	mkdir -p $@

Output/Deblur/%.table.qza Output/Deblur/%.seq.qza Output/Deblur/%.stat.qza: Output/Filter/%.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length 265 --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step05: Output/Deblur/whole.seq.qza
.PHONY += step05

# Step 06 (Taxonomy information)
Output/Taxonomy:
	mkdir -p $@

Output/Taxonomy/silva.qza: | Output/Taxonomy build.log
	wget "https://data.qiime2.org/2020.6/common/silva-138-99-nb-classifier.qza" -O $@

Output/Taxonomy/gg.qza: | Output/Taxonomy build.log
	wget "https://data.qiime2.org/2020.6/common/gg-13-8-99-nb-classifier.qza" -O $@

Output/Taxonomy/%.DADA2.silva.qza: Output/Taxonomy/silva.qza Output/DADA2/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.DADA2.gg.qza: Output/Taxonomy/gg.qza Output/DADA2/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.Deblur.silva.qza: Output/Taxonomy/silva.qza Output/Deblur/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.Deblur.gg.qza: Output/Taxonomy/gg.qza Output/Deblur/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step06: Output/Taxonomy/whole.DADA2.silva.qza Output/Taxonomy/whole.DADA2.gg.qza Output/Taxonomy/whole.Deblur.silva.qza Output/Taxonomy/whole.Deblur.gg.qza
.PHONY += step06

# Step 07 (Build Metadata)
Output/Metadata:
	mkdir -p $@

Output/Metadata/whole.tsv: Python/step07.py Data/Metadata/샘플정보.xlsx $(wildcard Data/RawData/*.fastq.gz) | Output/Metadata build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

step07: Output/Metadata/whole.tsv
.PHONY += step07

# Step 08 (Run ANCOM)
Output/ANCOM:
	mkdir -p $@

Output/ANCOM/%.DADA2.gg.qza: Output/DADA2/%.table.qza Output/Taxonomy/%.DADA2.gg.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.DADA2.silva.qza: Output/DADA2/%.table.qza Output/Taxonomy/%.DADA2.silva.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.Deblur.gg.qza: Output/Deblur/%.table.qza Output/Taxonomy/%.Deblur.gg.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.Deblur.silva.qza: Output/Deblur/%.table.qza Output/Taxonomy/%.Deblur.silva.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.CV.qza: Output/ANCOM/%.qza Output/Metadata/whole.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[site] IN ('Cervix', 'Vagina')" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.mouth.qza: Output/ANCOM/%.qza Output/Metadata/whole.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[site] IN ('Mouth')" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.pseudocount.qza: Output/ANCOM/%.qza | build.log
	$(DOCKER) qiime composition add-pseudocount --i-table $(addprefix /,$(word 1,$^)) --o-composition-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.site.qzv: Output/ANCOM/%.pseudocount.qza Output/Metadata/whole.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "site" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.site_premature.qzv: Output/ANCOM/%.pseudocount.qza Output/Metadata/whole.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "site_premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.CV.premature.qzv: Output/ANCOM/%.CV.pseudocount.qza Output/Metadata/whole.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.mouth.premature.qzv: Output/ANCOM/%.mouth.pseudocount.qza Output/Metadata/whole.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%/ancom.tsv: Output/ANCOM/%.qzv | build.log
	$(DOCKER) rm -rfv $(addprefix /,$(@D))
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D))

step08: Output/ANCOM/whole.DADA2.gg.site/ancom.tsv Output/ANCOM/whole.DADA2.silva.site/ancom.tsv Output/ANCOM/whole.Deblur.gg.site/ancom.tsv Output/ANCOM/whole.Deblur.silva.site/ancom.tsv
.PHONY += step08

# Step 09 (Get TSV)
Output/TSV:
	mkdir -p $@

Output/TSV/raw.DADA2.gg.tsv: Output/DADA2/whole.table.qza Output/Taxonomy/whole.DADA2.gg.qza | Output/TSV build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path /Output/TSV/$(RANDOM)
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 2,$^)) --output-path /Output/TSV/$(RANDOM)
	$(DOCKER) sed --in-place "1c#OTU ID\ttaxonomy\tconfidence" /Output/TSV/$(RANDOM)/taxonomy.tsv
	$(DOCKER) biom add-metadata --input-fp /Output/TSV/$(RANDOM)/feature-table.biom --observation-metadata-fp /Output/TSV/$(RANDOM)/taxonomy.tsv --output-fp /Output/TSV/$(RANDOM)/tmp.biom --sc-separated "taxonomy"
	$(DOCKER) biom convert --input-fp /Output/TSV/$(RANDOM)/tmp.biom --output-fp $(addprefix /,$@) --to-tsv --process-obs-metadata "taxonomy" --tsv-metadata-formatter "sc_separated" --header-key "taxonomy"
	$(DOCKER) rm -rfv /Output/TSV/$(RANDOM)

Output/TSV/raw.DADA2.silva.tsv: Output/DADA2/whole.table.qza Output/Taxonomy/whole.DADA2.silva.qza | Output/TSV build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path /Output/TSV/$(RANDOM)
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 2,$^)) --output-path /Output/TSV/$(RANDOM)
	$(DOCKER) sed --in-place "1c#OTU ID\ttaxonomy\tconfidence" /Output/TSV/$(RANDOM)/taxonomy.tsv
	$(DOCKER) biom add-metadata --input-fp /Output/TSV/$(RANDOM)/feature-table.biom --observation-metadata-fp /Output/TSV/$(RANDOM)/taxonomy.tsv --output-fp /Output/TSV/$(RANDOM)/tmp.biom --sc-separated "taxonomy"
	$(DOCKER) biom convert --input-fp /Output/TSV/$(RANDOM)/tmp.biom --output-fp $(addprefix /,$@) --to-tsv --process-obs-metadata "taxonomy" --tsv-metadata-formatter "sc_separated" --header-key "taxonomy"
	$(DOCKER) rm -rfv /Output/TSV/$(RANDOM)

Output/TSV/raw.Deblur.gg.tsv: Output/Deblur/whole.table.qza Output/Taxonomy/whole.Deblur.gg.qza | Output/TSV build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path /Output/TSV/$(RANDOM)
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 2,$^)) --output-path /Output/TSV/$(RANDOM)
	$(DOCKER) sed --in-place "1c#OTU ID\ttaxonomy\tconfidence" /Output/TSV/$(RANDOM)/taxonomy.tsv
	$(DOCKER) biom add-metadata --input-fp /Output/TSV/$(RANDOM)/feature-table.biom --observation-metadata-fp /Output/TSV/$(RANDOM)/taxonomy.tsv --output-fp /Output/TSV/$(RANDOM)/tmp.biom --sc-separated "taxonomy"
	$(DOCKER) biom convert --input-fp /Output/TSV/$(RANDOM)/tmp.biom --output-fp $(addprefix /,$@) --to-tsv --process-obs-metadata "taxonomy" --tsv-metadata-formatter "sc_separated" --header-key "taxonomy"
	$(DOCKER) rm -rfv /Output/TSV/$(RANDOM)

Output/TSV/raw.Deblur.silva.tsv: Output/Deblur/whole.table.qza Output/Taxonomy/whole.Deblur.silva.qza | Output/TSV build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path /Output/TSV/$(RANDOM)
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 2,$^)) --output-path /Output/TSV/$(RANDOM)
	$(DOCKER) sed --in-place "1c#OTU ID\ttaxonomy\tconfidence" /Output/TSV/$(RANDOM)/taxonomy.tsv
	$(DOCKER) biom add-metadata --input-fp /Output/TSV/$(RANDOM)/feature-table.biom --observation-metadata-fp /Output/TSV/$(RANDOM)/taxonomy.tsv --output-fp /Output/TSV/$(RANDOM)/tmp.biom --sc-separated "taxonomy"
	$(DOCKER) biom convert --input-fp /Output/TSV/$(RANDOM)/tmp.biom --output-fp $(addprefix /,$@) --to-tsv --process-obs-metadata "taxonomy" --tsv-metadata-formatter "sc_separated" --header-key "taxonomy"
	$(DOCKER) rm -rfv /Output/TSV/$(RANDOM)

step09: Output/TSV/raw.DADA2.gg.tsv Output/TSV/raw.DADA2.silva.tsv Output/TSV/raw.Deblur.gg.tsv Output/TSV/raw.Deblur.silva.tsv
.PHONY += step09

# Step 10 (Select & Clearify TSV)
Output/TSV/site.%.tar.gz: Python/step10.py Output/TSV/raw.%.tsv Output/ANCOM/whole.%.site/ancom.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^) $(addprefix /,$@)

Output/TSV/site_premature.%.tar.gz: Python/step10.py Output/TSV/raw.%.tsv Output/ANCOM/whole.%.site_premature/ancom.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^) $(addprefix /,$@)

Output/TSV/whole.%.tar.gz: Python/step10-2.py Output/TSV/raw.%.tsv
	$(DOCKER) python3 $(addprefix /,$^) $(addprefix /,$@)

step10: Output/TSV/site.DADA2.gg.tar.gz Output/TSV/site.DADA2.silva.tar.gz Output/TSV/site.Deblur.gg.tar.gz Output/TSV/site.Deblur.silva.tar.gz
.PHONY += step10

# Step 11 (Making t-SNE)
Output/TSNE:
	mkdir -p $@

Output/TSNE/%.tar.gz: Python/step11.py Output/TSV/%.tar.gz | Output/TSNE build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpu $(CPUS)

step11: Output/TSNE/site.DADA2.gg.tar.gz Output/TSNE/site.DADA2.silva.tar.gz Output/TSNE/site.Deblur.gg.tar.gz Output/TSNE/site.Deblur.silva.tar.gz
.PHONY += step11

# Step 12 (Drawing t-SNE)
Output/TSNE/%.png: Python/step12.py Output/TSNE/%.tar.gz Output/Metadata/whole.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@)

step12: Output/TSNE/site.DADA2.gg.png Output/TSNE/site_premature.DADA2.gg.png Output/TSNE/site.DADA2.silva.png Output/TSNE/site_premature.DADA2.silva.png Output/TSNE/site.Deblur.gg.png Output/TSNE/site_premature.Deblur.gg.png Output/TSNE/site.Deblur.silva.png Output/TSNE/site_premature.Deblur.silva.png
.PHONY += step12

# Step 13 (Build stability file)
MOTHUR_DOCKER = docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) --workdir $(addprefix /,$(@D)) $(IMAGE_NAME)
Output/Stability:
	mkdir -p $@

Output/Stability/whole.tsv: Python/step13.py $(wildcard Data/RawData/*.fastq.gz) | Output/Stability build.log
	$(DOCKER) python3 $(addprefix /,$^) > $@

step13: Output/Stability/whole.tsv
.PHONY += step13

# Step 14 (Make contigs)
Output/Contigs:
	mkdir -p $@

Output/Contigs/whole.trim.contigs.fasta Output/Contigs/whole.contigs.groups: Output/Stability/whole.tsv | Output/Contigs build.log
	$(MOTHUR_DOCKER) mothur "#set.dir(input=.); set.dir(output=.); make.contigs(file=$(addprefix /,$(word 1,$^)), processors=$(CPUS));" 1> $@.stdout 2> $@.stderr

step14: Output/Contigs/whole.trim.contigs.fasta
.PHONY += step14

# Step 15 (Screen sequence)
Output/Screen:
	mkdir -p $@

Output/Screen/whole.trim.contigs.good.fasta Output/Screen/whole.contigs.good.groups: Output/Contigs/whole.trim.contigs.fasta | Output/Screen build.log
	$(MOTHUR_DOCKER) mothur "#set.dir(input=.); set.dir(output=.); screen.seqs(fasta=$(addprefix /,$(word 1,$^)), processors=$(CPUS));"

step15: Output/Screen/whole.trim.contigs.good.fasta
.PHONY += step15

# Step 16 (Unique sequence)
Output/Unique:
	mkdir -p $@

Output/Unique/whole.trim.contigs.good.unique.fasta Output/Unique/whole.trim.contigs.good.names: Output/Screen/whole.trim.contigs.good.fasta | Output/Unique build.log
	$(MOTHUR_DOCKER) mothur "#set.dir(input=.); set.dir(output=.); unique.seqs(fasta=$(addprefix /,$(word 1,$^)));"

step16: Output/Unique/whole.trim.contigs.good.unique.fasta
.PHONY += step16

# Step 17 (Count sequence)
Output/Count:
	mkdir -p $@

Output/Count/whole.trim.contigs.good.count_table: Output/Unique/whole.trim.contigs.good.names | Output/Count build.log
	$(MOTHUR_DOCKER) mothur "#set.dir(input=.); set.dir(output=.); count.seqs(name=$(addprefix /,$(word 1,$^)));"

step17: Output/Count/whole.trim.contigs.good.count_table
.PHONY += step17

# Step 18 (Pseudo-PCR)
Output/PCR:
	mkdir -p $@

Output/PCR/silva.seed.tar.gz: | Output/PCR
	wget "https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138.tgz" -O $@

Output/PCR/silva.seed_v138.align Output/PCR/silva.seed_v138.tax: Output/PCR/silva.seed.tar.gz
	tar -zxvf $<
	touch $@

Output/PCR/silva.seed_v138.summary: Output/PCR/silva.seed_v138.align | build.log
	$(MOTHUR_DOCKER) mothur "#set.dir(input=.); set.dir(output=.); summary.seqs(fasta=$(addprefix /,$(word 1,$^)), processors=$(CPUS));"

Output/PCR/silva.seed_v138.scrap.pcr.align: Output/PCR/silva.seed_v138.align Output/PCR/silva.seed_v138.tax | build.log
	@echo "You have to check the length manually!!"
	$(MOTHUR_DOCKER) mothur "#set.dir(input=.); set.dir(output=.); pcr.seqs(fasta=$(addprefix /,$(word 1,$^)), taxonomy=$(addprefix /,$(word 2,$^)), start=1044, end=43116, keepdots=F, processors=$(CPUS));"

step18: Output/PCR/silva.seed_v138.scrap.pcr.align
.PHONY += step18

# Step 19 (Align to reference)
Output/Align:
	mkdir -p $@

Output/Align/whole.trim.contigs.good.unique.align: Output/Unique/whole.trim.contigs.good.unique.fasta Output/PCR/silva.seed_v138.scrap.pcr.align | Output/Align build.log
	$(MOTHUR_DOCKER) mothur "#set.dir(input=.); set.dir(output=.); align.seqs(fasta=$(addprefix /,$(word 1,$^)), reference=$(addprefix /,$(word 2,$^)), processors=$(CPUS));"

Output/Align/whole.trim.contigs.good.unique.summary: Output/Align/whole.trim.contigs.good.unique.align Output/Count/whole.trim.contigs.good.count_table | build.log
	$(MOTHUR_DOCKER) mothur "#set.dir(input=.); set.dir(output=.); summary.seqs(fasta=$(addprefix /,$(word 1,$^)), count=$(addprefix /,$(word 2,$^)), processors=$(CPUS));"

Output/Align/whole.trim.contigs.good.unique.good.align: Output/Align/whole.trim.contigs.good.unique.align Output/Count/whole.trim.contigs.good.count_table Output/Align/whole.trim.contigs.good.unique.summary | build.log
	@echo "You have to check the length manually!!"
	$(MOTHUR_DOCKER) mothur "#set.dir(input=.); set.dir(output=.); screen.seqs(fasta=$(addprefix /,$(word 1,$^)), count=$(addprefix /,$(word 2,$^)), summary=$(addprefix /,$(word 3,$^)), start=6388, end=25316, maxhomop=8, processors=$(CPUS));"

Output/Align/whole.trim.contigs.good.unique.good.filter.fasta: Output/Align/whole.trim.contigs.good.unique.good.align | build.log
	$(MOTHUR_DOCKER) mothur "#set.dir(input=.); set.dir(output=.); filter.seqs(fasta=$(addprefix /,$(word 1,$^)), vertical=T, trump=., processors=$(CPUS));"

Output/Align/whole.trim.contigs.good.unique.good.filter.unique.fasta Output/Align/whole.trim.contigs.good.unique.good.filter.count_table: Output/Align/whole.trim.contigs.good.unique.good.filter.fasta Output/Count/whole.trim.contigs.good.count_table
	$(MOTHUR_DOCKER) mothur "#set.dir(input=.); set.dir(output=.); unique.seqs(fasta=$(addprefix /,$(word 1,$^)), count=$(addprefix /,$(word 2,$^)))"

step19: Output/Align/whole.trim.contigs.good.unique.good.filter.unique.fasta
.PHONY += step19

# Step 20 (Chimera filtering)
Output/Chimera:
	mkdir -p $@

Output/Chimera/test: Output/Align/whole.trim.contigs.good.unique.good.filter.unique.fasta Output/Align/whole.trim.contigs.good.unique.good.filter.count_table | Output/Chimera build.log
	$(MOTHUR_DOCKER) mothur "#set.dir(input=.); set.dir(output=.); pre.cluster(fasta=$(addprefix /,$(word 1,$^)), count=$(addprefix /,$(word 2,$^)), diffs=2, processors=$(CPUS));"
# TODO: Here you are start

# Step 21 (Read & clearify raw TSV into pandas)
Output/TSV/ready.%.tar.gz: Python/step21.py Output/TSV/raw.%.tsv | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@)

step21: Output/TSV/ready.DADA2.gg.tar.gz Output/TSV/ready.DADA2.silva.tar.gz Output/TSV/ready.Deblur.gg.tar.gz Output/TSV/ready.Deblur.silva.tar.gz
.PHONY += step21

# Step 22 (Making & Draw t-SNE)
Output/Step22:
	mkdir -p $@

Output/Step22/%.tar.gz: Python/step11.py Output/TSV/%.tar.gz | Output/Step22 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) --cpu $(CPUS)

Output/Step22/%.png: Python/step22.py Output/Step22/%.tar.gz | build.log
	$(DOCKER) python3 $(addprefix /,$^ $@)

step22: Output/Step22/ready.DADA2.gg.png Output/Step22/ready.DADA2.silva.png Output/Step22/ready.Deblur.gg.png Output/Step22/ready.Deblur.silva.png
.PHONY += step22

# Step 23 (Draw histogram)
Output/Step23:
	mkdir -p $@

Output/Step23/age.png: Python/step23.py Data/Metadata/샘플정보.xlsx | Output/Step23 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) "나이"

Output/Step23/weight.png: Python/step23.py Data/Metadata/샘플정보.xlsx | Output/Step23 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) "임신전 몸무게"

Output/Step23/BMI.png: Python/step23.py Data/Metadata/샘플정보.xlsx | Output/Step23 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) "BMI"

step23: Output/Step23/age.png Output/Step23/weight.png Output/Step23/BMI.png
.PHONY += step23

# Step 24 (Draw histogram with parsing)
Output/Step24:
	mkdir -p $@

Output/Step24/weight_now.png: Python/step24.py Data/Metadata/샘플정보.xlsx | Output/Step24 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) "체중/키" 0 "/"

Output/Step24/height.png: Python/step24.py Data/Metadata/샘플정보.xlsx | Output/Step24 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) "체중/키" 1 "/"

Output/Step24/weeks.png: Python/step24.py Data/Metadata/샘플정보.xlsx | Output/Step24 build.log
	$(DOCKER) python3 $(addprefix /,$^ $@) "분만주수" 0 "+"

step24: Output/Step24/weight_now.png Output/Step24/height.png Output/Step24/weeks.png
.PHONY += step24
