# Variables
DATE := $(shell date "+%y%m%d")
RANDOM := $(shell bash -c 'echo $$RANDOM')
IMAGE_NAME = premature:latest
SHELL := /bin/bash
CPUS = 10
MEMS = 100G
PWD := $(shell pwd)
TOOLS =
DOCKER = docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) $(IMAGE_NAME)
PYTHON = /opt/conda/envs/qiime2-2022.2/bin/python3 -B
BASH = /bin/bash
R = /opt/conda/envs/qiime2-2022.2/bin/Rscript --vanilla

# Options
VOLUME_OPTS = --volume $(abspath Output):/Output --volume $(abspath Data):/Data:ro --volume $(abspath Python):/Python:ro --volume $(abspath Bash):/Bash:ro --volume $(abspath R):/R:ro --tmpfs /tmpfs
RUN_OPTS = --tty --cpus="$(CPUS)" --memory="$(MEMS)"

# General
all:
.PHONY: all

log Output Docker/Tools:
	mkdir $@

# Docker
TOOLS += $(wildcard Docker/*)
build.log: Docker/Dockerfile $(TOOLS) | log Output
	rm -fv $@
	docker images | grep $(IMAGE_NAME) && docker rmi $(IMAGE_NAME) || true
	docker build --rm --tag $(IMAGE_NAME) $(<D) | tee $@

build: build.log
.PHONY: build

interactive: build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) --interactive $(IMAGE_NAME) $(BASH) || true
.PHONY: interactive

delete: build.log
	docker rmi $(IMAGE_NAME)
	rm -fv build.log
.PHONY: delete

stop:
	docker rm $(CONTAINER_NAME)

# SGE
run: | log Output
	echo -e "#!$(BASH)\nmake -C $(PWD) latest" > tmp.sh
	sbatch --chdir=$(abspath .) --cpus-per-task=$(CPUS) --error="$(abspath log)/%x.e%A" --job-name="Premature_$(DATE)" --mail-type=ALL --mail-user="230@fumire.moe" --mem=$(MEMS) --output="$(abspath log)/%x.o%A" --export=ALL tmp.sh
.PHONY: run

# Actual
latest: step01 step02 step03 step04 step05 step06 step07 step08 step09 step21 step22 step26 step28 step29 step30 step31 step32 step33 step35 step36 step37 step38 step39 step40 step41 step42 step43 step44 step45 step46 step47 step48 step49 step50 step51 step52 step53 step54 step55 step56 step57 step58 step59 step60 step61 step62 step63 step64 step65 step66 step67 step68 step69 step70 step71
.PHONY: latest

# Step 01 (Make manifest files)
Output/Manifest:
	mkdir -p $@

Output/Manifest/first.tsv: Python/step01.py $(wildcard Data/RawData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^) > $@

Output/Manifest/second.tsv: Python/step01.py $(wildcard Data/HN00145920/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^) > $@

Output/Manifest/EBI.tsv: Python/step01-2.py $(wildcard Data/EBIData/ena_files/*/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^) > $@

Output/Manifest/HMP.tsv: Python/step01-3.py $(wildcard Data/HMPData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^) > $@

Output/Manifest/Stool.tsv: Python/step01-4.py Data/Metadata/readable.xlsx $(wildcard Data/StoolData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^) $(filter %.xlsx,$^)) $(foreach f,$(wildcard Data/StoolData/*.fastq.gz),--stool $(addprefix /,$f)) > $@

Output/Manifest/everything.tsv: Python/step01-5.py Data/Metadata/readable.xlsx $(wildcard Data/RawData/*.fastq.gz) $(wildcard Data/HN00145920/*.fastq.gz) $(wildcard Data/HN00148220/*.fastq.gz) $(wildcard Data/StoolData/*.fastq.gz) | Output/Manifest build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^) $(filter %.xlsx,$^)) $(foreach f,$(wildcard Data/RawData/*.fastq.gz),--first $(addprefix /,$f)) $(foreach f,$(wildcard Data/HN00145920/*.fastq.gz),--second $(addprefix /,$f)) $(foreach f,$(wildcard Data/HN00148220/*.fastq.gz),--third $(addprefix /,$f)) > $@

step01: Output/Manifest/Stool.tsv Output/Manifest/everything.tsv
.PHONY: step01

# Step 02 (De-multiplexing)
Output/Demux:
	mkdir -p $@

Output/Demux/first.qza: Output/Manifest/first.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/second.qza: Output/Manifest/second.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/EBI.qza: Output/Manifest/EBI.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[SequencesWithQuality]" --input-format "SingleEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/HMP.qza: Output/Manifest/HMP.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/Stool.qza: Output/Manifest/Stool.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/everything.qza: Output/Manifest/everything.tsv | Output/Demux build.log
	$(DOCKER) qiime tools import --type "SampleData[PairedEndSequencesWithQuality]" --input-format "PairedEndFastqManifestPhred33V2" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Demux/%.qzv: Output/Demux/%.qza | build.log
	$(DOCKER) qiime demux summarize --i-data $(addprefix /,$<) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step02: Output/Demux/Stool.qza Output/Demux/Stool.qzv Output/Demux/everything.qza Output/Demux/everything.qzv
.PHONY: step02

# Step 03 (Quality filter)
Output/Filter:
	mkdir -p $@

Output/Filter/%.seq.qza Output/Filter/%.stat.qza: Output/Demux/%.qza | Output/Filter build.log
	$(DOCKER) qiime quality-filter q-score --i-demux $(addprefix /,$<) --o-filtered-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-filter-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step03: Output/Filter/Stool.seq.qza Output/Filter/everything.seq.qza
.PHONY: step03

# Step 04 (DADA2 Denoising)
Output/DADA2:
	mkdir -p $@

Output/DADA2/first.table.qza Output/DADA2/first.seq.qza Output/DADA2/first.stat.qza: Output/Demux/first.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "300" --p-trunc-len-r "265" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/second.table.qza Output/DADA2/second.seq.qza Output/DADA2/second.stat.qza: Output/Demux/second.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "300" --p-trunc-len-r "222" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/EBI.table.qza Output/DADA2/EBI.seq.qza Output/DADA2/EBI.stat.qza: Output/Demux/EBI.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-single --i-demultiplexed-seqs $(addprefix /,$<) --p-trunc-len "150" --p-n-threads $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/HMP.table.qza Output/DADA2/HMP.seq.qza Output/DADA2/HMP.stat.qza: Output/Demux/HMP.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "278" --p-trunc-len-r "226" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/Stool.table.qza Output/DADA2/Stool.seq.qza Output/DADA2/Stool.stat.qza: Output/Demux/Stool.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "250" --p-trunc-len-r "251" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/DADA2/everything.table.qza Output/DADA2/everything.seq.qza Output/DADA2/everything.stat.qza: Output/Demux/everything.qza | Output/DADA2 build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f "300" --p-trunc-len-r "245" --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step04: Output/DADA2/Stool.seq.qza Output/DADA2/everything.seq.qza
.PHONY: step04

# Step 05 (Deblur Denoising)
Output/Deblur:
	mkdir -p $@

Output/Deblur/first.table.qza Output/Deblur/first.seq.qza Output/Deblur/first.stat.qza: Output/Filter/first.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "265" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/second.table.qza Output/Deblur/second.seq.qza Output/Deblur/second.stat.qza: Output/Filter/second.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "222" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/EBI.table.qza Output/Deblur/EBI.seq.qza Output/Deblur/EBI.stat.qza: Output/Filter/EBI.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "150" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/HMP.table.qza Output/Deblur/HMP.seq.qza Output/Deblur/HMP.stat.qza: Output/Filter/HMP.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "226" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/Stool.table.qza Output/Deblur/Stool.seq.qza Output/Deblur/Stool.stat.qza: Output/Filter/Stool.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "250" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

Output/Deblur/everything.table.qza Output/Deblur/everything.seq.qza Output/Deblur/everything.stat.qza: Output/Filter/everything.seq.qza | Output/Deblur build.log
	@echo "You have to check the length manually!!"
	$(DOCKER) qiime deblur denoise-16S --i-demultiplexed-seqs $(addprefix /,$<) --p-trim-length "245" --p-jobs-to-start $(CPUS) --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step05: Output/Deblur/Stool.seq.qza Output/Deblur/everything.seq.qza
.PHONY: step05

# Step 06 (Taxonomy information)
Output/Taxonomy:
	mkdir -p $@

Output/Taxonomy/silva.qza: | Output/Taxonomy build.log
	wget "https://data.qiime2.org/2021.4/common/silva-138-99-nb-classifier.qza" -O $@
	touch $@

Output/Taxonomy/gg.qza: | Output/Taxonomy build.log
	wget "https://data.qiime2.org/2021.4/common/gg-13-8-99-nb-classifier.qza" -O $@
	touch $@

Output/Taxonomy/homd.fasta: | Output/Taxonomy
	wget "http://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/V15.22/HOMD_16S_rRNA_RefSeq_V15.22.fasta" -O $@
	touch $@

Output/Taxonomy/homd.txt: | Output/Taxonomy
	wget "http://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/V15.22/HOMD_16S_rRNA_RefSeq_V15.22.qiime.taxonomy" -O $@
	touch $@

Output/Taxonomy/homd_fasta.qza: Output/Taxonomy/homd.fasta | build.log
	$(DOCKER) qiime tools import --type "FeatureData[Sequence]" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/homd_taxo.qza: Output/Taxonomy/homd.txt | build.log
	$(DOCKER) qiime tools import --type "FeatureData[Taxonomy]" --input-format "HeaderlessTSVTaxonomyFormat" --input-path $(addprefix /,$<) --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/homd_seqs.qza: Output/Taxonomy/homd_fasta.qza | build.log
	$(DOCKER) qiime feature-classifier extract-reads --i-sequences $(addprefix /,$<) --p-f-primer "GTGCCAGCMGCCGCGGTAA" --p-r-primer "GGACTACHVGGGTWTCTAAT" --o-reads $(addprefix /,$@) --p-n-jobs $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/homd.qza: Output/Taxonomy/homd_seqs.qza Output/Taxonomy/homd_taxo.qza | build.log
	$(DOCKER) qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads $(addprefix /,$(word 1,$^)) --i-reference-taxonomy $(addprefix /,$(word 2,$^)) --o-classifier $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.DADA2.silva.qza: Output/Taxonomy/silva.qza Output/DADA2/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.DADA2.gg.qza: Output/Taxonomy/gg.qza Output/DADA2/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.DADA2.homd.qza: Output/Taxonomy/homd.qza Output/DADA2/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.Deblur.silva.qza: Output/Taxonomy/silva.qza Output/Deblur/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.Deblur.gg.qza: Output/Taxonomy/gg.qza Output/Deblur/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Taxonomy/%.Deblur.homd.qza: Output/Taxonomy/homd.qza Output/Deblur/%.seq.qza | build.log
	$(DOCKER) qiime feature-classifier classify-sklearn --i-classifier $(addprefix /,$(word 1,$^)) --i-reads $(addprefix /,$(word 2,$^)) --p-n-jobs $(CPUS) --o-classification $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step06-reference: Output/Taxonomy/silva.qza Output/Taxonomy/gg.qza Output/Taxonomy/homd.fasta Output/Taxonomy/homd.txt Output/Taxonomy/homd_fasta.qza Output/Taxonomy/homd_taxo.qza Output/Taxonomy/homd_seqs.qza Output/Taxonomy/homd.qza
step06-mapping: Output/Taxonomy/everything.DADA2.gg.qza Output/Taxonomy/everything.DADA2.silva.qza Output/Taxonomy/everything.DADA2.homd.qza Output/Taxonomy/everything.Deblur.gg.qza Output/Taxonomy/everything.Deblur.silva.qza Output/Taxonomy/everything.Deblur.homd.qza Output/Taxonomy/Stool.DADA2.gg.qza Output/Taxonomy/Stool.DADA2.silva.qza Output/Taxonomy/Stool.DADA2.homd.qza Output/Taxonomy/Stool.Deblur.gg.qza Output/Taxonomy/Stool.Deblur.silva.qza Output/Taxonomy/Stool.Deblur.homd.qza Output/Taxonomy/EBI.DADA2.gg.qza Output/Taxonomy/EBI.DADA2.silva.qza Output/Taxonomy/EBI.DADA2.homd.qza Output/Taxonomy/EBI.Deblur.gg.qza Output/Taxonomy/EBI.Deblur.silva.qza Output/Taxonomy/EBI.Deblur.homd.qza Output/Taxonomy/HMP.DADA2.gg.qza Output/Taxonomy/HMP.DADA2.silva.qza Output/Taxonomy/HMP.DADA2.homd.qza Output/Taxonomy/HMP.Deblur.gg.qza Output/Taxonomy/HMP.Deblur.silva.qza Output/Taxonomy/HMP.Deblur.homd.qza
step06: step06-reference step06-mapping
.PHONY: step06 step06-reference step06-mapping

# Step 07 (Build Metadata)
Output/Metadata:
	mkdir -p $@

Output/Metadata/everything.tsv: Python/step07-4.py Output/Manifest/everything.tsv Data/Metadata/metadata.xlsx | Output/Metadata build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^) > $@

Output/Metadata/Stool.tsv: Python/step07-4.py Output/Manifest/Stool.tsv Data/Metadata/metadata.xlsx | Output/Metadata build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^) > $@

step07: Output/Metadata/everything.tsv Output/Metadata/Stool.tsv
.PHONY: step07

# Step 08 (Run ANCOM)
Output/ANCOM:
	mkdir -p $@

Output/ANCOM/%.DADA2.gg.qza: Output/DADA2/%.table.qza Output/Taxonomy/%.DADA2.gg.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.DADA2.silva.qza: Output/DADA2/%.table.qza Output/Taxonomy/%.DADA2.silva.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.DADA2.homd.qza: Output/DADA2/%.table.qza Output/Taxonomy/%.DADA2.homd.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.Deblur.gg.qza: Output/Deblur/%.table.qza Output/Taxonomy/%.Deblur.gg.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.Deblur.silva.qza: Output/Deblur/%.table.qza Output/Taxonomy/%.Deblur.silva.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.Deblur.homd.qza: Output/Deblur/%.table.qza Output/Taxonomy/%.Deblur.homd.qza | Output/ANCOM build.log
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Mouth.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]='Mouth'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Cervix.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]='Cervix'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Vagina.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]='Vagina'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Neonate.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site] IN ('Neonate-1day', 'Neonate-3day', 'Neonate-5day')" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Neonate-1day.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]=='Neonate-1day'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Neonate-3day.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]=='Neonate-3day'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Neonate-5day.qza: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]=='Neonate-5day'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.pseudocount.qza: Output/ANCOM/%.qza | build.log
	$(DOCKER) qiime composition add-pseudocount --i-table $(addprefix /,$(word 1,$^)) --o-composition-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.PTB.qzv: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "Premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Detail_PTB.qzv: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "Detail Premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/everything.%.Simple_PTB.qzv: Output/ANCOM/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "Simple Premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%/ancom.tsv Output/ANCOM/%/data.tsv: Output/ANCOM/%.qzv | build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D))
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D)) 1> $(@D).stdout 2> $(@D).stderr

Output/ANCOM/%.pdf: Python/step08.py Output/ANCOM/%/ancom.tsv Output/ANCOM/%/data.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step08-DADA2-gg: Output/ANCOM/everything.DADA2.gg.qza Output/ANCOM/everything.DADA2.gg.pseudocount.qza Output/ANCOM/everything.DADA2.gg.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.gg.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.gg.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.gg.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.gg.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.gg.Mouth.qza Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.gg.Cervix.qza Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.gg.Vagina.qza Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate-1day.qza Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Neonate-1day.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate-3day.qza Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Neonate-3day.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate-3day.qza Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.qza Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.gg.Neonate-5day.pseudocount.Simple_PTB.pdf
step08-DADA2-silva: Output/ANCOM/everything.DADA2.silva.qza Output/ANCOM/everything.DADA2.silva.pseudocount.qza Output/ANCOM/everything.DADA2.silva.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.silva.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.silva.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.silva.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.silva.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.silva.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.silva.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.silva.Mouth.qza Output/ANCOM/everything.DADA2.silva.Mouth.pseudocount.qza Output/ANCOM/everything.DADA2.silva.Mouth.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.silva.Mouth.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Mouth.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.silva.Mouth.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.silva.Mouth.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Mouth.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.silva.Mouth.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.silva.Mouth.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Mouth.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.silva.Cervix.qza Output/ANCOM/everything.DADA2.silva.Cervix.pseudocount.qza Output/ANCOM/everything.DADA2.silva.Cervix.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.silva.Cervix.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Cervix.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.silva.Cervix.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.silva.Cervix.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Cervix.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.silva.Cervix.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.silva.Cervix.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Cervix.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.silva.Vagina.qza Output/ANCOM/everything.DADA2.silva.Vagina.pseudocount.qza Output/ANCOM/everything.DADA2.silva.Vagina.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.silva.Vagina.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Vagina.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.silva.Vagina.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.silva.Vagina.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Vagina.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.silva.Vagina.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.silva.Vagina.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Vagina.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.silva.Neonate-1day.qza Output/ANCOM/everything.DADA2.silva.Neonate-1day.pseudocount.qza Output/ANCOM/everything.DADA2.silva.Neonate-1day.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.silva.Neonate-1day.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Neonate-1day.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.silva.Neonate-1day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.silva.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Neonate-1day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.silva.Neonate-1day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.silva.Neonate-1day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Neonate-1day.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.silva.Neonate-3day.qza Output/ANCOM/everything.DADA2.silva.Neonate-3day.pseudocount.qza Output/ANCOM/everything.DADA2.silva.Neonate-3day.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.silva.Neonate-3day.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Neonate-3day.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.silva.Neonate-3day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.silva.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Neonate-3day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.silva.Neonate-3day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.silva.Neonate-3day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Neonate-3day.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.silva.Neonate-3day.qza Output/ANCOM/everything.DADA2.silva.Neonate-5day.pseudocount.qza Output/ANCOM/everything.DADA2.silva.Neonate-5day.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.silva.Neonate-5day.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Neonate-5day.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.silva.Neonate-5day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.silva.Neonate-5day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Neonate-5day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.silva.Neonate-5day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.silva.Neonate-5day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.silva.Neonate-5day.pseudocount.Simple_PTB.pdf
step08-DADA2-homd: Output/ANCOM/everything.DADA2.homd.qza Output/ANCOM/everything.DADA2.homd.pseudocount.qza Output/ANCOM/everything.DADA2.homd.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.homd.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.homd.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.homd.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.homd.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.homd.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.homd.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.homd.Mouth.qza Output/ANCOM/everything.DADA2.homd.Mouth.pseudocount.qza Output/ANCOM/everything.DADA2.homd.Mouth.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.homd.Mouth.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Mouth.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.homd.Mouth.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.homd.Mouth.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Mouth.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.homd.Mouth.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.homd.Mouth.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Mouth.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.homd.Cervix.qza Output/ANCOM/everything.DADA2.homd.Cervix.pseudocount.qza Output/ANCOM/everything.DADA2.homd.Cervix.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.homd.Cervix.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Cervix.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.homd.Cervix.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.homd.Cervix.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Cervix.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.homd.Cervix.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.homd.Cervix.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Cervix.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.homd.Vagina.qza Output/ANCOM/everything.DADA2.homd.Vagina.pseudocount.qza Output/ANCOM/everything.DADA2.homd.Vagina.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.homd.Vagina.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Vagina.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.homd.Vagina.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.homd.Vagina.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Vagina.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.homd.Vagina.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.homd.Vagina.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Vagina.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.homd.Neonate-1day.qza Output/ANCOM/everything.DADA2.homd.Neonate-1day.pseudocount.qza Output/ANCOM/everything.DADA2.homd.Neonate-1day.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.homd.Neonate-1day.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Neonate-1day.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.homd.Neonate-1day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.homd.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Neonate-1day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.homd.Neonate-1day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.homd.Neonate-1day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Neonate-1day.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.homd.Neonate-3day.qza Output/ANCOM/everything.DADA2.homd.Neonate-3day.pseudocount.qza Output/ANCOM/everything.DADA2.homd.Neonate-3day.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.homd.Neonate-3day.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Neonate-3day.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.homd.Neonate-3day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.homd.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Neonate-3day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.homd.Neonate-3day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.homd.Neonate-3day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Neonate-3day.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.DADA2.homd.Neonate-3day.qza Output/ANCOM/everything.DADA2.homd.Neonate-5day.pseudocount.qza Output/ANCOM/everything.DADA2.homd.Neonate-5day.pseudocount.PTB.qzv Output/ANCOM/everything.DADA2.homd.Neonate-5day.pseudocount.PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Neonate-5day.pseudocount.PTB.pdf Output/ANCOM/everything.DADA2.homd.Neonate-5day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.DADA2.homd.Neonate-5day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Neonate-5day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.DADA2.homd.Neonate-5day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.DADA2.homd.Neonate-5day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.DADA2.homd.Neonate-5day.pseudocount.Simple_PTB.pdf
step08-Deblur-gg: Output/ANCOM/everything.Deblur.gg.qza Output/ANCOM/everything.Deblur.gg.pseudocount.qza Output/ANCOM/everything.Deblur.gg.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.gg.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.gg.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.gg.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.gg.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.gg.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.gg.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.gg.Mouth.qza Output/ANCOM/everything.Deblur.gg.Mouth.pseudocount.qza Output/ANCOM/everything.Deblur.gg.Mouth.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.gg.Mouth.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Mouth.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.gg.Mouth.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.gg.Mouth.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Mouth.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.gg.Mouth.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.gg.Mouth.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Mouth.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.gg.Cervix.qza Output/ANCOM/everything.Deblur.gg.Cervix.pseudocount.qza Output/ANCOM/everything.Deblur.gg.Cervix.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.gg.Cervix.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Cervix.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.gg.Cervix.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.gg.Cervix.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Cervix.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.gg.Cervix.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.gg.Cervix.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Cervix.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.gg.Vagina.qza Output/ANCOM/everything.Deblur.gg.Vagina.pseudocount.qza Output/ANCOM/everything.Deblur.gg.Vagina.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.gg.Vagina.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Vagina.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.gg.Vagina.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.gg.Vagina.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Vagina.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.gg.Vagina.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.gg.Vagina.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Vagina.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.gg.Neonate-1day.qza Output/ANCOM/everything.Deblur.gg.Neonate-1day.pseudocount.qza Output/ANCOM/everything.Deblur.gg.Neonate-1day.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.gg.Neonate-1day.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Neonate-1day.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.gg.Neonate-1day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.gg.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Neonate-1day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.gg.Neonate-1day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.gg.Neonate-1day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Neonate-1day.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.gg.Neonate-3day.qza Output/ANCOM/everything.Deblur.gg.Neonate-3day.pseudocount.qza Output/ANCOM/everything.Deblur.gg.Neonate-3day.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.gg.Neonate-3day.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Neonate-3day.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.gg.Neonate-3day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.gg.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Neonate-3day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.gg.Neonate-3day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.gg.Neonate-3day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Neonate-3day.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.gg.Neonate-5day.qza Output/ANCOM/everything.Deblur.gg.Neonate-5day.pseudocount.qza Output/ANCOM/everything.Deblur.gg.Neonate-5day.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.gg.Neonate-5day.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Neonate-5day.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.gg.Neonate-5day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.gg.Neonate-5day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Neonate-5day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.gg.Neonate-5day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.gg.Neonate-5day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.gg.Neonate-5day.pseudocount.Simple_PTB.pdf
step08-Deblur-silva: Output/ANCOM/everything.Deblur.silva.qza Output/ANCOM/everything.Deblur.silva.pseudocount.qza Output/ANCOM/everything.Deblur.silva.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.silva.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.silva.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.silva.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.silva.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.silva.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.silva.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.silva.Mouth.qza Output/ANCOM/everything.Deblur.silva.Mouth.pseudocount.qza Output/ANCOM/everything.Deblur.silva.Mouth.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.silva.Mouth.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Mouth.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.silva.Mouth.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.silva.Mouth.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Mouth.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.silva.Mouth.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.silva.Mouth.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Mouth.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.silva.Cervix.qza Output/ANCOM/everything.Deblur.silva.Cervix.pseudocount.qza Output/ANCOM/everything.Deblur.silva.Cervix.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.silva.Cervix.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Cervix.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.silva.Cervix.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.silva.Cervix.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Cervix.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.silva.Cervix.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.silva.Cervix.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Cervix.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.silva.Vagina.qza Output/ANCOM/everything.Deblur.silva.Vagina.pseudocount.qza Output/ANCOM/everything.Deblur.silva.Vagina.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.silva.Vagina.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Vagina.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.silva.Vagina.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.silva.Vagina.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Vagina.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.silva.Vagina.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.silva.Vagina.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Vagina.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.silva.Neonate-1day.qza Output/ANCOM/everything.Deblur.silva.Neonate-1day.pseudocount.qza Output/ANCOM/everything.Deblur.silva.Neonate-1day.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.silva.Neonate-1day.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Neonate-1day.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.silva.Neonate-1day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.silva.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Neonate-1day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.silva.Neonate-1day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.silva.Neonate-1day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Neonate-1day.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.silva.Neonate-3day.qza Output/ANCOM/everything.Deblur.silva.Neonate-3day.pseudocount.qza Output/ANCOM/everything.Deblur.silva.Neonate-3day.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.silva.Neonate-3day.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Neonate-3day.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.silva.Neonate-3day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.silva.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Neonate-3day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.silva.Neonate-3day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.silva.Neonate-3day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Neonate-3day.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.silva.Neonate-3day.qza Output/ANCOM/everything.Deblur.silva.Neonate-5day.pseudocount.qza Output/ANCOM/everything.Deblur.silva.Neonate-5day.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.silva.Neonate-5day.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Neonate-5day.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.silva.Neonate-5day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.silva.Neonate-5day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Neonate-5day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.silva.Neonate-5day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.silva.Neonate-5day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.silva.Neonate-5day.pseudocount.Simple_PTB.pdf
step08-Deblur-homd: Output/ANCOM/everything.Deblur.homd.qza Output/ANCOM/everything.Deblur.homd.pseudocount.qza Output/ANCOM/everything.Deblur.homd.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.homd.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.homd.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.homd.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.homd.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.homd.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.homd.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.homd.Mouth.qza Output/ANCOM/everything.Deblur.homd.Mouth.pseudocount.qza Output/ANCOM/everything.Deblur.homd.Mouth.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.homd.Mouth.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Mouth.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.homd.Mouth.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.homd.Mouth.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Mouth.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.homd.Mouth.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.homd.Mouth.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Mouth.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.homd.Cervix.qza Output/ANCOM/everything.Deblur.homd.Cervix.pseudocount.qza Output/ANCOM/everything.Deblur.homd.Cervix.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.homd.Cervix.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Cervix.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.homd.Cervix.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.homd.Cervix.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Cervix.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.homd.Cervix.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.homd.Cervix.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Cervix.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.homd.Vagina.qza Output/ANCOM/everything.Deblur.homd.Vagina.pseudocount.qza Output/ANCOM/everything.Deblur.homd.Vagina.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.homd.Vagina.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Vagina.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.homd.Vagina.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.homd.Vagina.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Vagina.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.homd.Vagina.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.homd.Vagina.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Vagina.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.homd.Neonate-1day.qza Output/ANCOM/everything.Deblur.homd.Neonate-1day.pseudocount.qza Output/ANCOM/everything.Deblur.homd.Neonate-1day.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.homd.Neonate-1day.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Neonate-1day.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.homd.Neonate-1day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.homd.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Neonate-1day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.homd.Neonate-1day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.homd.Neonate-1day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Neonate-1day.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.homd.Neonate-3day.qza Output/ANCOM/everything.Deblur.homd.Neonate-3day.pseudocount.qza Output/ANCOM/everything.Deblur.homd.Neonate-3day.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.homd.Neonate-3day.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Neonate-3day.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.homd.Neonate-3day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.homd.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Neonate-3day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.homd.Neonate-3day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.homd.Neonate-3day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Neonate-3day.pseudocount.Simple_PTB.pdf Output/ANCOM/everything.Deblur.homd.Neonate-3day.qza Output/ANCOM/everything.Deblur.homd.Neonate-5day.pseudocount.qza Output/ANCOM/everything.Deblur.homd.Neonate-5day.pseudocount.PTB.qzv Output/ANCOM/everything.Deblur.homd.Neonate-5day.pseudocount.PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Neonate-5day.pseudocount.PTB.pdf Output/ANCOM/everything.Deblur.homd.Neonate-5day.pseudocount.Detail_PTB.qzv Output/ANCOM/everything.Deblur.homd.Neonate-5day.pseudocount.Detail_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Neonate-5day.pseudocount.Detail_PTB.pdf Output/ANCOM/everything.Deblur.homd.Neonate-5day.pseudocount.Simple_PTB.qzv Output/ANCOM/everything.Deblur.homd.Neonate-5day.pseudocount.Simple_PTB/data.tsv Output/ANCOM/everything.Deblur.homd.Neonate-5day.pseudocount.Simple_PTB.pdf
step08: step08-DADA2-gg step08-DADA2-silva step08-DADA2-homd step08-Deblur-gg step08-Deblur-silva step08-Deblur-homd
.PHONY: step08-DADA2-gg step08-DADA2-silva step08-DADA2-homd step08-Deblur-gg step08-Deblur-silva step08-Deblur-homd step08

# Step 09 (Get TSV)
Output/TSV:
	mkdir -p $@

Output/TSV/everything.DADA2.%.uncorrected.tsv: Bash/export_TSV.bash Output/DADA2/everything.table.qza Output/Taxonomy/everything.DADA2.%.qza | Output/TSV build.log
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/everything.Deblur.%.uncorrected.tsv: Bash/export_TSV.bash Output/Deblur/everything.table.qza Output/Taxonomy/everything.Deblur.%.qza | Output/TSV build.log
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/Macrogen.uncorrected.tsv: Python/step09.py Data/Macrogen/First/otu_table_shared_with_tax_assignment.xlsx Data/Macrogen/Second+Third/Mouth/otu_table_shared_with_tax_assignment.xlsx Data/Macrogen/Second+Third/Vagina/otu_table_shared_with_tax_assignment.xlsx Data/Macrogen/Second+Third/Neonate-1day/otu_table_shared_with_tax_assignment.xlsx Data/Macrogen/Second+Third/Neonate-3day/otu_table_shared_with_tax_assignment.xlsx Data/Macrogen/Second+Third/Neonate-5day/otu_table_shared_with_tax_assignment.xlsx Data/Metadata/readable.xlsx | Output/TSV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/%.log.tsv: Python/step09-1.py Output/TSV/%.uncorrected.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/TSV/%.tsv: R/remove_batch.R Output/TSV/%.uncorrected.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step09-uncorrected: Output/TSV/everything.DADA2.gg.uncorrected.tsv Output/TSV/everything.DADA2.silva.uncorrected.tsv Output/TSV/everything.DADA2.homd.uncorrected.tsv Output/TSV/everything.Deblur.gg.uncorrected.tsv Output/TSV/everything.Deblur.silva.uncorrected.tsv Output/TSV/everything.Deblur.homd.uncorrected.tsv Output/TSV/Macrogen.uncorrected.tsv
step09-log: Output/TSV/everything.DADA2.gg.log.tsv Output/TSV/everything.DADA2.silva.log.tsv Output/TSV/everything.DADA2.homd.log.tsv Output/TSV/everything.Deblur.gg.log.tsv Output/TSV/everything.Deblur.silva.log.tsv Output/TSV/everything.Deblur.homd.log.tsv Output/TSV/Macrogen.tsv
step09-limma: Output/TSV/everything.DADA2.gg.tsv Output/TSV/everything.DADA2.silva.tsv Output/TSV/everything.DADA2.homd.tsv Output/TSV/everything.Deblur.gg.tsv Output/TSV/everything.Deblur.silva.tsv Output/TSV/everything.Deblur.homd.tsv Output/TSV/Macrogen.tsv
step09: step09-uncorrected step09-log step09-limma
.PHONY: step09-uncorrected step09-log step09-limma step09

# Step 21 (Read & clearify first TSV into pandas)
Output/TSV_Proportion:
	mkdir -p $@

Output/TSV_Proportion/%.tsv: Python/step21.py Output/TSV/%.tsv | build.log Output/TSV_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step21-uncorrected: Output/TSV_Proportion/everything.DADA2.gg.uncorrected.tsv Output/TSV_Proportion/everything.DADA2.silva.uncorrected.tsv Output/TSV_Proportion/everything.DADA2.homd.uncorrected.tsv Output/TSV_Proportion/everything.Deblur.gg.uncorrected.tsv Output/TSV_Proportion/everything.Deblur.silva.uncorrected.tsv Output/TSV_Proportion/everything.Deblur.homd.uncorrected.tsv
step21-log: Output/TSV_Proportion/everything.DADA2.gg.log.tsv Output/TSV_Proportion/everything.DADA2.silva.log.tsv Output/TSV_Proportion/everything.DADA2.homd.log.tsv Output/TSV_Proportion/everything.Deblur.gg.log.tsv Output/TSV_Proportion/everything.Deblur.silva.log.tsv Output/TSV_Proportion/everything.Deblur.homd.log.tsv
step21-corrected: Output/TSV_Proportion/everything.DADA2.gg.tsv Output/TSV_Proportion/everything.DADA2.silva.tsv Output/TSV_Proportion/everything.DADA2.homd.tsv Output/TSV_Proportion/everything.Deblur.gg.tsv Output/TSV_Proportion/everything.Deblur.silva.tsv Output/TSV_Proportion/everything.Deblur.homd.tsv
step21: step21-uncorrected step21-log step21-corrected
.PHONY: step21-uncorrected step21-log step21-corrected step21

# Step 22 (Making & Draw t-SNE)
Output/Step22 Output/Step22_Proportion:
	mkdir -p $@

Output/Step22/everything.%.tsv: Python/step11.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step22 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step22_Proportion/everything.%.tsv: Python/step11.py Output/TSV_Proportion/everything.%.tsv Output/Metadata/everything.tsv | Output/Step22_Proportion build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step22/%.tar: Python/step22-3.py Output/Step22/%.tsv | Output/Step22 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step22_Proportion/%.tar: Python/step22-3.py Output/Step22_Proportion/%.tsv | Output/Step22 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step22-uncorrected: Output/Step22/everything.DADA2.gg.uncorrected.tsv Output/Step22/everything.DADA2.gg.uncorrected.tar Output/Step22/everything.DADA2.silva.uncorrected.tsv Output/Step22/everything.DADA2.silva.uncorrected.tar Output/Step22/everything.DADA2.homd.uncorrected.tsv Output/Step22/everything.DADA2.homd.uncorrected.tar Output/Step22/everything.Deblur.gg.uncorrected.tsv Output/Step22/everything.Deblur.gg.uncorrected.tar Output/Step22/everything.Deblur.silva.uncorrected.tsv Output/Step22/everything.Deblur.silva.uncorrected.tar Output/Step22/everything.Deblur.homd.uncorrected.tsv Output/Step22/everything.Deblur.homd.uncorrected.tar
step22-log: Output/Step22/everything.DADA2.gg.log.tsv Output/Step22/everything.DADA2.gg.log.tar Output/Step22/everything.DADA2.silva.log.tsv Output/Step22/everything.DADA2.silva.log.tar Output/Step22/everything.DADA2.homd.log.tsv Output/Step22/everything.DADA2.homd.log.tar Output/Step22/everything.Deblur.gg.log.tsv Output/Step22/everything.Deblur.gg.log.tar Output/Step22/everything.Deblur.silva.log.tsv Output/Step22/everything.Deblur.silva.log.tar Output/Step22/everything.Deblur.homd.log.tsv Output/Step22/everything.Deblur.homd.log.tar
step22-corrected: Output/Step22/everything.DADA2.gg.tsv Output/Step22/everything.DADA2.gg.tar Output/Step22/everything.DADA2.silva.tsv Output/Step22/everything.DADA2.silva.tar Output/Step22/everything.DADA2.homd.tsv Output/Step22/everything.DADA2.homd.tar Output/Step22/everything.Deblur.gg.tsv Output/Step22/everything.Deblur.gg.tar Output/Step22/everything.Deblur.silva.tsv Output/Step22/everything.Deblur.silva.tar Output/Step22/everything.Deblur.homd.tsv Output/Step22/everything.Deblur.homd.tar
step22-uncorrected-proportion: Output/Step22_Proportion/everything.DADA2.gg.uncorrected.tsv Output/Step22_Proportion/everything.DADA2.gg.uncorrected.tar Output/Step22_Proportion/everything.DADA2.silva.uncorrected.tsv Output/Step22_Proportion/everything.DADA2.silva.uncorrected.tar Output/Step22_Proportion/everything.DADA2.homd.uncorrected.tsv Output/Step22_Proportion/everything.DADA2.homd.uncorrected.tar Output/Step22_Proportion/everything.Deblur.gg.uncorrected.tsv Output/Step22_Proportion/everything.Deblur.gg.uncorrected.tar Output/Step22_Proportion/everything.Deblur.silva.uncorrected.tsv Output/Step22_Proportion/everything.Deblur.silva.uncorrected.tar Output/Step22_Proportion/everything.Deblur.homd.uncorrected.tsv Output/Step22_Proportion/everything.Deblur.homd.uncorrected.tar
step22-log-proportion: Output/Step22_Proportion/everything.DADA2.gg.log.tsv Output/Step22_Proportion/everything.DADA2.gg.log.tar Output/Step22_Proportion/everything.DADA2.silva.log.tsv Output/Step22_Proportion/everything.DADA2.silva.log.tar Output/Step22_Proportion/everything.DADA2.homd.log.tsv Output/Step22_Proportion/everything.DADA2.homd.log.tar Output/Step22_Proportion/everything.Deblur.gg.log.tsv Output/Step22_Proportion/everything.Deblur.gg.log.tar Output/Step22_Proportion/everything.Deblur.silva.log.tsv Output/Step22_Proportion/everything.Deblur.silva.log.tar Output/Step22_Proportion/everything.Deblur.homd.log.tsv Output/Step22_Proportion/everything.Deblur.homd.log.tar
step22-corrected-proportion: Output/Step22_Proportion/everything.DADA2.gg.tsv Output/Step22_Proportion/everything.DADA2.gg.tar Output/Step22_Proportion/everything.DADA2.silva.tsv Output/Step22_Proportion/everything.DADA2.silva.tar Output/Step22_Proportion/everything.DADA2.homd.tsv Output/Step22_Proportion/everything.DADA2.homd.tar Output/Step22_Proportion/everything.Deblur.gg.tsv Output/Step22_Proportion/everything.Deblur.gg.tar Output/Step22_Proportion/everything.Deblur.silva.tsv Output/Step22_Proportion/everything.Deblur.silva.tar Output/Step22_Proportion/everything.Deblur.homd.tsv Output/Step22_Proportion/everything.Deblur.homd.tar
step22: step22-uncorrected step22-log step22-corrected step22-uncorrected-proportion step22-log-proportion step22-corrected-proportion
.PHONY: step22-uncorrected step22-log step22-corrected step22-uncorrected-proportion step22-log-proportion step22-corrected-proportion step22

# Step 26 (RandomForest Classifier)
Output/Step26 Output/Step26_Proportion:
	mkdir -p $@

Output/Step26/RF.%.tar: Python/step26.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step26 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step26/RF-two.%.tar: Python/step26-1.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step26 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step26_Proportion/RF.%.tar: Python/step26.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step26_Proportion build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step26_Proportion/RF-two.%.tar: Python/step26-1.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step26_Proportion build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step26-uncorrected: Output/Step26/RF.DADA2.gg.uncorrected.tar Output/Step26/RF.DADA2.silva.uncorrected.tar Output/Step26/RF.DADA2.homd.uncorrected.tar Output/Step26/RF.Deblur.gg.uncorrected.tar Output/Step26/RF.Deblur.silva.uncorrected.tar Output/Step26/RF.Deblur.homd.uncorrected.tar Output/Step26/RF-two.DADA2.gg.uncorrected.tar Output/Step26/RF-two.DADA2.silva.uncorrected.tar Output/Step26/RF-two.DADA2.homd.uncorrected.tar Output/Step26/RF-two.Deblur.gg.uncorrected.tar Output/Step26/RF-two.Deblur.silva.uncorrected.tar Output/Step26/RF-two.Deblur.homd.uncorrected.tar
step26-corrected: Output/Step26/RF.DADA2.gg.tar Output/Step26/RF.DADA2.silva.tar Output/Step26/RF.DADA2.homd.tar Output/Step26/RF.Deblur.gg.tar Output/Step26/RF.Deblur.silva.tar Output/Step26/RF.Deblur.homd.tar Output/Step26/RF-two.DADA2.gg.tar Output/Step26/RF-two.DADA2.silva.tar Output/Step26/RF-two.DADA2.homd.tar Output/Step26/RF-two.Deblur.gg.tar Output/Step26/RF-two.Deblur.silva.tar Output/Step26/RF-two.Deblur.homd.tar
step26-uncorrected-Proportion: Output/Step26_Proportion/RF.DADA2.gg.uncorrected.tar Output/Step26_Proportion/RF.DADA2.silva.uncorrected.tar Output/Step26_Proportion/RF.DADA2.homd.uncorrected.tar Output/Step26_Proportion/RF.Deblur.gg.uncorrected.tar Output/Step26_Proportion/RF.Deblur.silva.uncorrected.tar Output/Step26_Proportion/RF.Deblur.homd.uncorrected.tar Output/Step26_Proportion/RF-two.DADA2.gg.uncorrected.tar Output/Step26_Proportion/RF-two.DADA2.silva.uncorrected.tar Output/Step26_Proportion/RF-two.DADA2.homd.uncorrected.tar Output/Step26_Proportion/RF-two.Deblur.gg.uncorrected.tar Output/Step26_Proportion/RF-two.Deblur.silva.uncorrected.tar Output/Step26_Proportion/RF-two.Deblur.homd.uncorrected.tar
step26-corrected-Proportion: Output/Step26_Proportion/RF.DADA2.gg.tar Output/Step26_Proportion/RF.DADA2.silva.tar Output/Step26_Proportion/RF.DADA2.homd.tar Output/Step26_Proportion/RF.Deblur.gg.tar Output/Step26_Proportion/RF.Deblur.silva.tar Output/Step26_Proportion/RF.Deblur.homd.tar Output/Step26_Proportion/RF-two.DADA2.gg.tar Output/Step26_Proportion/RF-two.DADA2.silva.tar Output/Step26_Proportion/RF-two.DADA2.homd.tar Output/Step26_Proportion/RF-two.Deblur.gg.tar Output/Step26_Proportion/RF-two.Deblur.silva.tar Output/Step26_Proportion/RF-two.Deblur.homd.tar
step26: step26-uncorrected step26-corrected step26-uncorrected-Proportion step26-corrected-Proportion
.PHONY: step26-uncorrected step26-corrected step26-uncorrected-Proportion step26-corrected-Proportion step26

# Step 27 (DecisionTree Classifier)
Output/Step27:
	mkdir -p $@

Output/Step27/DT.%.tar: Python/step27.py Output/TSV/everything.%.tar.gz Output/TSV/EBI.%.tar.gz Output/TSV/HMP.%.tar.gz Output/Metadata/everything.tsv | Output/Step27 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpu $(CPUS) 1> $@.stdout 2> $@.stderr

step27: Output/Step27/DT.DADA2.gg.tar Output/Step27/DT.DADA2.silva.tar Output/Step27/DT.Deblur.gg.tar Output/Step27/DT.Deblur.silva.tar
.PHONY: step27

# Step 28 (Phylogenetic diversity analysis)
Output/Step28:
	mkdir -p $@

Output/Step28/%.DADA2.aligned-rep-seqs.qza Output/Step28/%.DADA2.masked-aligned-rep-seqs.qza Output/Step28/%.DADA2.unrooted-tree.qza Output/Step28/%.DADA2.rooted-tree.qza: Output/DADA2/%.seq.qza | Output/Step28 build.log
	$(DOCKER) qiime phylogeny align-to-tree-mafft-fasttree --i-sequences $(addprefix /,$<) --o-alignment $(addprefix /,$(basename $(basename $@))).aligned-rep-seqs.qza --o-masked-alignment $(addprefix /,$(basename $(basename $@))).masked-aligned-rep-seqs.qza --o-tree $(addprefix /,$(basename $(basename $@))).unrooted-tree.qza --o-rooted-tree $(addprefix /,$(basename $(basename $@))).rooted-tree.qza --p-n-threads $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step28/%.Deblur.aligned-rep-seqs.qza Output/Step28/%.Deblur.masked-aligned-rep-seqs.qza Output/Step28/%.Deblur.unrooted-tree.qza Output/Step28/%.Deblur.rooted-tree.qza: Output/Deblur/%.seq.qza | Output/Step28 build.log
	$(DOCKER) qiime phylogeny align-to-tree-mafft-fasttree --i-sequences $(addprefix /,$<) --o-alignment $(addprefix /,$(basename $(basename $@))).aligned-rep-seqs.qza --o-masked-alignment $(addprefix /,$(basename $(basename $@))).masked-aligned-rep-seqs.qza --o-tree $(addprefix /,$(basename $(basename $@))).unrooted-tree.qza --o-rooted-tree $(addprefix /,$(basename $(basename $@))).rooted-tree.qza --p-n-threads $(CPUS) 1> $@.stdout 2> $@.stderr

step28: Output/Step28/everything.DADA2.rooted-tree.qza
.PHONY: step28

# Step 29 (Core phylogenetic metrics)
Output/Step29 Output/Step29/DADA2 Output/Step29/Deblur:
	mkdir -p $@

Output/Step29/%.DADA2.table.qzv: Output/DADA2/%.table.qza Output/Metadata/%.tsv | Output/Step29 build.log
	$(DOCKER) qiime feature-table summarize --i-table $(addprefix /,$(word 1,$^)) --m-sample-metadata-file $(addprefix /,$(word 2,$^)) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step29/%.DADA2.seq.qzv: Output/DADA2/%.seq.qza | Output/Step29 build.log
	$(DOCKER) qiime feature-table tabulate-seqs --i-data $(addprefix /,$<) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step29/%.Deblur.table.qzv: Output/Deblur/%.table.qza Output/Metadata/%.tsv | Output/Step29 build.log
	$(DOCKER) qiime feature-table summarize --i-table $(addprefix /,$(word 1,$^)) --m-sample-metadata-file $(addprefix /,$(word 2,$^)) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step29/%.Deblur.seq.qzv: Output/Deblur/%.seq.qza | Output/Step29 build.log
	$(DOCKER) qiime feature-table tabulate-seqs --i-data $(addprefix /,$<) --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step29/DADA2/everything.faith_pd_vector.qza Output/Step29/DADA2/everything.unweighted_unifrac_distance_matrix.qza Output/Step29/DADA2/everything.bray_curtis_pcoa_results.qza Output/Step29/DADA2/everything.shannon_vector.qza Output/Step29/DADA2/everything.rarefied_table.qza Output/Step29/DADA2/everything.weighted_unifrac_distance_matrix.qza Output/Step29/DADA2/everything.jaccard_pcoa_results.qza Output/Step29/DADA2/everything.weighted_unifrac_pcoa_results.qza Output/Step29/DADA2/everything.observed_features_vector.qza Output/Step29/DADA2/everything.jaccard_distance_matrix.qza Output/Step29/DADA2/everything.evenness_vector.qza Output/Step29/DADA2/everything.bray_curtis_distance_matrix.qza Output/Step29/DADA2/everything.unweighted_unifrac_pcoa_results.qza Output/Step29/DADA2/everything.unweighted_unifrac_emperor.qzv Output/Step29/DADA2/everything.jaccard_emperor.qzv Output/Step29/DADA2/everything.bray_curtis_emperor.qzv Output/Step29/DADA2/everything.weighted_unifrac_emperor.qzv: Output/Step28/everything.DADA2.rooted-tree.qza Output/DADA2/everything.table.qza Output/Metadata/everything.tsv | Output/Step29/DADA2 build.log
	@echo "You have to check the sequencing depth manually!!"
	$(DOCKER) qiime diversity core-metrics-phylogenetic --i-phylogeny $(addprefix /,$(word 1,$^)) --i-table $(addprefix /,$(word 2,$^)) --m-metadata-file $(addprefix /,$(word 3,$^)) --p-sampling-depth "1046" --o-rarefied-table $(addprefix /,$(basename $(basename $@))).rarefied_table.qza --o-faith-pd-vector $(addprefix /,$(basename $(basename $@))).faith_pd_vector.qza --o-observed-features-vector $(addprefix /,$(basename $(basename $@))).observed_features_vector.qza --o-shannon-vector $(addprefix /,$(basename $(basename $@))).shannon_vector.qza --o-evenness-vector $(addprefix /,$(basename $(basename $@))).evenness_vector.qza --o-unweighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_distance_matrix.qza --o-weighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).weighted_unifrac_distance_matrix.qza --o-jaccard-distance-matrix $(addprefix /,$(basename $(basename $@))).jaccard_distance_matrix.qza --o-bray-curtis-distance-matrix $(addprefix /,$(basename $(basename $@))).bray_curtis_distance_matrix.qza --o-unweighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_pcoa_results.qza --o-weighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).weighted_unifrac_pcoa_results.qza --o-unweighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_emperor.qzv --o-weighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).weighted_unifrac_emperor.qzv --o-jaccard-emperor $(addprefix /,$(basename $(basename $@))).jaccard_emperor.qzv --o-bray-curtis-emperor $(addprefix /,$(basename $(basename $@))).bray_curtis_emperor.qzv --o-jaccard-pcoa-results $(addprefix /,$(basename $(basename $@))).jaccard_pcoa_results.qza --o-bray-curtis-pcoa-results $(addprefix /,$(basename $(basename $@))).bray_curtis_pcoa_results.qza --p-n-jobs-or-threads $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step29/DADA2/Stool.faith_pd_vector.qza Output/Step29/DADA2/Stool.unweighted_unifrac_distance_matrix.qza Output/Step29/DADA2/Stool.bray_curtis_pcoa_results.qza Output/Step29/DADA2/Stool.shannon_vector.qza Output/Step29/DADA2/Stool.rarefied_table.qza Output/Step29/DADA2/Stool.weighted_unifrac_distance_matrix.qza Output/Step29/DADA2/Stool.jaccard_pcoa_results.qza Output/Step29/DADA2/Stool.weighted_unifrac_pcoa_results.qza Output/Step29/DADA2/Stool.observed_features_vector.qza Output/Step29/DADA2/Stool.jaccard_distance_matrix.qza Output/Step29/DADA2/Stool.evenness_vector.qza Output/Step29/DADA2/Stool.bray_curtis_distance_matrix.qza Output/Step29/DADA2/Stool.unweighted_unifrac_pcoa_results.qza Output/Step29/DADA2/Stool.unweighted_unifrac_emperor.qzv Output/Step29/DADA2/Stool.jaccard_emperor.qzv Output/Step29/DADA2/Stool.bray_curtis_emperor.qzv Output/Step29/DADA2/Stool.weighted_unifrac_emperor.qzv: Output/Step28/Stool.DADA2.rooted-tree.qza Output/DADA2/Stool.table.qza Output/Metadata/Stool.tsv | Output/Step29/DADA2 build.log
	@echo "You have to check the sequencing depth manually!!"
	$(DOCKER) qiime diversity core-metrics-phylogenetic --i-phylogeny $(addprefix /,$(word 1,$^)) --i-table $(addprefix /,$(word 2,$^)) --m-metadata-file $(addprefix /,$(word 3,$^)) --p-sampling-depth "7992" --o-rarefied-table $(addprefix /,$(basename $(basename $@))).rarefied_table.qza --o-faith-pd-vector $(addprefix /,$(basename $(basename $@))).faith_pd_vector.qza --o-observed-features-vector $(addprefix /,$(basename $(basename $@))).observed_features_vector.qza --o-shannon-vector $(addprefix /,$(basename $(basename $@))).shannon_vector.qza --o-evenness-vector $(addprefix /,$(basename $(basename $@))).evenness_vector.qza --o-unweighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_distance_matrix.qza --o-weighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).weighted_unifrac_distance_matrix.qza --o-jaccard-distance-matrix $(addprefix /,$(basename $(basename $@))).jaccard_distance_matrix.qza --o-bray-curtis-distance-matrix $(addprefix /,$(basename $(basename $@))).bray_curtis_distance_matrix.qza --o-unweighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_pcoa_results.qza --o-weighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).weighted_unifrac_pcoa_results.qza --o-unweighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_emperor.qzv --o-weighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).weighted_unifrac_emperor.qzv --o-jaccard-emperor $(addprefix /,$(basename $(basename $@))).jaccard_emperor.qzv --o-bray-curtis-emperor $(addprefix /,$(basename $(basename $@))).bray_curtis_emperor.qzv --o-jaccard-pcoa-results $(addprefix /,$(basename $(basename $@))).jaccard_pcoa_results.qza --o-bray-curtis-pcoa-results $(addprefix /,$(basename $(basename $@))).bray_curtis_pcoa_results.qza --p-n-jobs-or-threads $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step29/Deblur/everything.faith_pd_vector.qza Output/Step29/Deblur/everything.unweighted_unifrac_distance_matrix.qza Output/Step29/Deblur/everything.bray_curtis_pcoa_results.qza Output/Step29/Deblur/everything.shannon_vector.qza Output/Step29/Deblur/everything.rarefied_table.qza Output/Step29/Deblur/everything.weighted_unifrac_distance_matrix.qza Output/Step29/Deblur/everything.jaccard_pcoa_results.qza Output/Step29/Deblur/everything.weighted_unifrac_pcoa_results.qza Output/Step29/Deblur/everything.observed_features_vector.qza Output/Step29/Deblur/everything.jaccard_distance_matrix.qza Output/Step29/Deblur/everything.evenness_vector.qza Output/Step29/Deblur/everything.bray_curtis_distance_matrix.qza Output/Step29/Deblur/everything.unweighted_unifrac_pcoa_results.qza Output/Step29/Deblur/everything.unweighted_unifrac_emperor.qzv Output/Step29/Deblur/everything.jaccard_emperor.qzv Output/Step29/Deblur/everything.bray_curtis_emperor.qzv Output/Step29/Deblur/everything.weighted_unifrac_emperor.qzv: Output/Step28/everything.Deblur.rooted-tree.qza Output/Deblur/everything.table.qza Output/Metadata/everything.tsv | Output/Step29/Deblur build.log
	@echo "You have to check the sequencing depth manually!!"
	$(DOCKER) qiime diversity core-metrics-phylogenetic --i-phylogeny $(addprefix /,$(word 1,$^)) --i-table $(addprefix /,$(word 2,$^)) --m-metadata-file $(addprefix /,$(word 3,$^)) --p-sampling-depth "4864" --o-rarefied-table $(addprefix /,$(basename $(basename $@))).rarefied_table.qza --o-faith-pd-vector $(addprefix /,$(basename $(basename $@))).faith_pd_vector.qza --o-observed-features-vector $(addprefix /,$(basename $(basename $@))).observed_features_vector.qza --o-shannon-vector $(addprefix /,$(basename $(basename $@))).shannon_vector.qza --o-evenness-vector $(addprefix /,$(basename $(basename $@))).evenness_vector.qza --o-unweighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_distance_matrix.qza --o-weighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).weighted_unifrac_distance_matrix.qza --o-jaccard-distance-matrix $(addprefix /,$(basename $(basename $@))).jaccard_distance_matrix.qza --o-bray-curtis-distance-matrix $(addprefix /,$(basename $(basename $@))).bray_curtis_distance_matrix.qza --o-unweighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_pcoa_results.qza --o-weighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).weighted_unifrac_pcoa_results.qza --o-unweighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_emperor.qzv --o-weighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).weighted_unifrac_emperor.qzv --o-jaccard-emperor $(addprefix /,$(basename $(basename $@))).jaccard_emperor.qzv --o-bray-curtis-emperor $(addprefix /,$(basename $(basename $@))).bray_curtis_emperor.qzv --o-jaccard-pcoa-results $(addprefix /,$(basename $(basename $@))).jaccard_pcoa_results.qza --o-bray-curtis-pcoa-results $(addprefix /,$(basename $(basename $@))).bray_curtis_pcoa_results.qza --p-n-jobs-or-threads $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step29/Deblur/Stool.faith_pd_vector.qza Output/Step29/Deblur/Stool.unweighted_unifrac_distance_matrix.qza Output/Step29/Deblur/Stool.bray_curtis_pcoa_results.qza Output/Step29/Deblur/Stool.shannon_vector.qza Output/Step29/Deblur/Stool.rarefied_table.qza Output/Step29/Deblur/Stool.weighted_unifrac_distance_matrix.qza Output/Step29/Deblur/Stool.jaccard_pcoa_results.qza Output/Step29/Deblur/Stool.weighted_unifrac_pcoa_results.qza Output/Step29/Deblur/Stool.observed_features_vector.qza Output/Step29/Deblur/Stool.jaccard_distance_matrix.qza Output/Step29/Deblur/Stool.evenness_vector.qza Output/Step29/Deblur/Stool.bray_curtis_distance_matrix.qza Output/Step29/Deblur/Stool.unweighted_unifrac_pcoa_results.qza Output/Step29/Deblur/Stool.unweighted_unifrac_emperor.qzv Output/Step29/Deblur/Stool.jaccard_emperor.qzv Output/Step29/Deblur/Stool.bray_curtis_emperor.qzv Output/Step29/Deblur/Stool.weighted_unifrac_emperor.qzv: Output/Step28/Stool.Deblur.rooted-tree.qza Output/Deblur/Stool.table.qza Output/Metadata/Stool.tsv | Output/Step29/Deblur build.log
	@echo "You have to check the sequencing depth manually!!"
	$(DOCKER) qiime diversity core-metrics-phylogenetic --i-phylogeny $(addprefix /,$(word 1,$^)) --i-table $(addprefix /,$(word 2,$^)) --m-metadata-file $(addprefix /,$(word 3,$^)) --p-sampling-depth "4645" --o-rarefied-table $(addprefix /,$(basename $(basename $@))).rarefied_table.qza --o-faith-pd-vector $(addprefix /,$(basename $(basename $@))).faith_pd_vector.qza --o-observed-features-vector $(addprefix /,$(basename $(basename $@))).observed_features_vector.qza --o-shannon-vector $(addprefix /,$(basename $(basename $@))).shannon_vector.qza --o-evenness-vector $(addprefix /,$(basename $(basename $@))).evenness_vector.qza --o-unweighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_distance_matrix.qza --o-weighted-unifrac-distance-matrix $(addprefix /,$(basename $(basename $@))).weighted_unifrac_distance_matrix.qza --o-jaccard-distance-matrix $(addprefix /,$(basename $(basename $@))).jaccard_distance_matrix.qza --o-bray-curtis-distance-matrix $(addprefix /,$(basename $(basename $@))).bray_curtis_distance_matrix.qza --o-unweighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_pcoa_results.qza --o-weighted-unifrac-pcoa-results $(addprefix /,$(basename $(basename $@))).weighted_unifrac_pcoa_results.qza --o-unweighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).unweighted_unifrac_emperor.qzv --o-weighted-unifrac-emperor $(addprefix /,$(basename $(basename $@))).weighted_unifrac_emperor.qzv --o-jaccard-emperor $(addprefix /,$(basename $(basename $@))).jaccard_emperor.qzv --o-bray-curtis-emperor $(addprefix /,$(basename $(basename $@))).bray_curtis_emperor.qzv --o-jaccard-pcoa-results $(addprefix /,$(basename $(basename $@))).jaccard_pcoa_results.qza --o-bray-curtis-pcoa-results $(addprefix /,$(basename $(basename $@))).bray_curtis_pcoa_results.qza --p-n-jobs-or-threads $(CPUS) 1> $@.stdout 2> $@.stderr

step29-depth: Output/Step29/everything.DADA2.table.qzv Output/Step29/everything.DADA2.seq.qzv Output/Step29/everything.Deblur.table.qzv
#  Output/Step29/everything.Deblur.seq.qzv Output/Step29/Stool.DADA2.table.qzv Output/Step29/Stool.DADA2.seq.qzv Output/Step29/Stool.Deblur.table.qzv Output/Step29/Stool.Deblur.seq.qzv
step29: Output/Step29/DADA2/everything.faith_pd_vector.qza Output/Step29/Deblur/everything.faith_pd_vector.qza
# Output/Step29/DADA2/Stool.faith_pd_vector.qza Output/Step29/Deblur/Stool.faith_pd_vector.qza
.PHONY: step29-depth step29

# Step 30 (Alpha-diversity)
Output/Step30 Output/Step30/DADA2 Output/Step30/Deblur:
	mkdir -p $@

Output/Step30/DADA2/%.tsv: Output/Step29/DADA2/%_vector.qza | Output/Step30/DADA2 build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D)) 1> $@.stdout 2> $@.stderr
	$(DOCKER) mv -f $(addprefix /,$(@D))/alpha-diversity.tsv $(addprefix /,$@)

Output/Step30/Deblur/%.tsv: Output/Step29/Deblur/%_vector.qza | Output/Step30/Deblur build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D)) 1> $@.stdout 2> $@.stderr
	$(DOCKER) mv -f $(addprefix /,$(@D))/alpha-diversity.tsv $(addprefix /,$@)

Output/Step30/everything.DADA2.tar: Python/step30.py Output/Step30/DADA2/everything.evenness.tsv Output/Step30/DADA2/everything.faith_pd.tsv Output/Step30/DADA2/everything.observed_features.tsv Output/Step30/DADA2/everything.shannon.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step30/First.DADA2.tar: Python/step30.py Output/Step30/DADA2/everything.evenness.tsv Output/Step30/DADA2/everything.faith_pd.tsv Output/Step30/DADA2/everything.observed_features.tsv Output/Step30/DADA2/everything.shannon.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --first 1> $@.stdout 2> $@.stderr

Output/Step30/Second_Third.DADA2.tar: Python/step30.py Output/Step30/DADA2/everything.evenness.tsv Output/Step30/DADA2/everything.faith_pd.tsv Output/Step30/DADA2/everything.observed_features.tsv Output/Step30/DADA2/everything.shannon.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --second 1> $@.stdout 2> $@.stderr

Output/Step30/everything.Deblur.tar: Python/step30.py Output/Step30/Deblur/everything.evenness.tsv Output/Step30/Deblur/everything.faith_pd.tsv Output/Step30/Deblur/everything.observed_features.tsv Output/Step30/Deblur/everything.shannon.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step30/First.Deblur.tar: Python/step30.py Output/Step30/Deblur/everything.evenness.tsv Output/Step30/Deblur/everything.faith_pd.tsv Output/Step30/Deblur/everything.observed_features.tsv Output/Step30/Deblur/everything.shannon.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --first 1> $@.stdout 2> $@.stderr

Output/Step30/Second_Third.Deblur.tar: Python/step30.py Output/Step30/Deblur/everything.evenness.tsv Output/Step30/Deblur/everything.faith_pd.tsv Output/Step30/Deblur/everything.observed_features.tsv Output/Step30/Deblur/everything.shannon.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --second 1> $@.stdout 2> $@.stderr

Output/Step30/Macrogen.tar: Python/step30-1.py Output/TSV/Macrogen.tsv Output/Metadata/everything.tsv | Output/Step30 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step30/First.Macrogen.tar: Python/step30-1.py Output/TSV/Macrogen.tsv Output/Metadata/everything.tsv | Output/Step30 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --first 1> $@.stdout 2> $@.stderr

Output/Step30/Second_Third.Macrogen.tar: Python/step30-1.py Output/TSV/Macrogen.tsv Output/Metadata/everything.tsv | Output/Step30 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --second 1> $@.stdout 2> $@.stderr

step30-DADA2: Output/Step30/DADA2/everything.evenness.tsv Output/Step30/DADA2/everything.faith_pd.tsv Output/Step30/DADA2/everything.observed_features.tsv Output/Step30/DADA2/everything.shannon.tsv Output/Step30/everything.DADA2.tar Output/Step30/First.DADA2.tar Output/Step30/Second_Third.DADA2.tar
step30-Deblur: Output/Step30/Deblur/everything.evenness.tsv Output/Step30/Deblur/everything.faith_pd.tsv Output/Step30/Deblur/everything.observed_features.tsv Output/Step30/Deblur/everything.shannon.tsv Output/Step30/everything.Deblur.tar Output/Step30/First.Deblur.tar Output/Step30/Second_Third.Deblur.tar
step30-Macrogen: Output/Step30/Macrogen.tar Output/Step30/First.Macrogen.tar Output/Step30/Second_Third.Macrogen.tar
step30: step30-DADA2 step30-Deblur step30-Macrogen
.PHONY: step30-DADA2 step30-Deblur step30-Macrogen step30

# Step 31 (Violin plots from ANCOM)
Output/ANCOM/everything.DADA2.gg.%.tar: Python/step31.py Output/ANCOM/everything.DADA2.gg.%/ancom.tsv Output/TSV/everything.DADA2.gg.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step31: Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB.tar Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.PTB.tar Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.PTB.tar Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.PTB.tar Output/ANCOM/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB.tar Output/ANCOM/everything.DADA2.gg.Cervix.pseudocount.Detail_PTB.tar Output/ANCOM/everything.DADA2.gg.Vagina.pseudocount.Detail_PTB.tar Output/ANCOM/everything.DADA2.gg.Baby.pseudocount.Detail_PTB.tar
.PHONY: step31

# Step 32 (Beta-diversity)
Output/Step32 Output/Step32/Macrogen:
	mkdir -p $@

Output/Step32/everything.%.tar: Python/step32.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step32
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step32/Macrogen/everything.%.tar: Python/step32-1.py Output/TSV/Macrogen.tsv Output/Metadata/everything.tsv | Output/Step32/Macrogen build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --beta "$*" 1> $@.stdout 2> $@.stderr

step32-uncorrected: Output/Step32/everything.DADA2.gg.uncorrected.tar Output/Step32/everything.DADA2.silva.uncorrected.tar Output/Step32/everything.DADA2.homd.uncorrected.tar Output/Step32/everything.Deblur.gg.uncorrected.tar Output/Step32/everything.Deblur.silva.uncorrected.tar Output/Step32/everything.Deblur.homd.uncorrected.tar
step32-corrected: Output/Step32/everything.DADA2.gg.tar Output/Step32/everything.DADA2.silva.tar Output/Step32/everything.DADA2.homd.tar Output/Step32/everything.Deblur.gg.tar Output/Step32/everything.Deblur.silva.tar Output/Step32/everything.Deblur.homd.tar
step32-Macrogen: Output/Step32/Macrogen/everything.unweighted_unifrac.tar Output/Step32/Macrogen/First.unweighted_unifrac.tar Output/Step32/Macrogen/Second_Third.unweighted_unifrac.tar Output/Step32/Macrogen/everything.weighted_unifrac.tar Output/Step32/Macrogen/First.weighted_unifrac.tar Output/Step32/Macrogen/Second_Third.weighted_unifrac.tar
step32: step32-uncorrected step32-corrected
.PHONY: step32-uncorrected step32-corrected step32-Macrogen step32

# Step 33 (LEfSe analysis)
Output/Step33 Output/Step33/Default Output/Step33/Cladogram:
	mkdir -p $@
	chmod 777 $@

Output/Step33/%.DADA2.gg.Cervix.tsv Output/Step33/%.DADA2.gg.Mouth.tsv Output/Step33/%.DADA2.gg.Neonate-3day.tsv Output/Step33/%.DADA2.gg.Neonate-5day.tsv Output/Step33/%.DADA2.gg.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.gg.uncorrected.tsv Output/Metadata/%.tsv | Output/Step33 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/%.DADA2.silva.Cervix.tsv Output/Step33/%.DADA2.silva.Mouth.tsv Output/Step33/%.DADA2.silva.Neonate-3day.tsv Output/Step33/%.DADA2.silva.Neonate-5day.tsv Output/Step33/%.DADA2.silva.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.silva.uncorrected.tsv Output/Metadata/%.tsv | Output/Step33 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/%.DADA2.homd.Cervix.tsv Output/Step33/%.DADA2.homd.Mouth.tsv Output/Step33/%.DADA2.homd.Neonate-3day.tsv Output/Step33/%.DADA2.homd.Neonate-5day.tsv Output/Step33/%.DADA2.homd.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.homd.uncorrected.tsv Output/Metadata/%.tsv | Output/Step33 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/%.Deblur.gg.Cervix.tsv Output/Step33/%.Deblur.gg.Mouth.tsv Output/Step33/%.Deblur.gg.Neonate-3day.tsv Output/Step33/%.Deblur.gg.Neonate-5day.tsv Output/Step33/%.Deblur.gg.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.gg.uncorrected.tsv Output/Metadata/%.tsv | Output/Step33 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/%.Deblur.silva.Cervix.tsv Output/Step33/%.Deblur.silva.Mouth.tsv Output/Step33/%.Deblur.silva.Neonate-3day.tsv Output/Step33/%.Deblur.silva.Neonate-5day.tsv Output/Step33/%.Deblur.silva.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.silva.uncorrected.tsv Output/Metadata/%.tsv | Output/Step33 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/%.Deblur.homd.Cervix.tsv Output/Step33/%.Deblur.homd.Mouth.tsv Output/Step33/%.Deblur.homd.Neonate-3day.tsv Output/Step33/%.Deblur.homd.Neonate-5day.tsv Output/Step33/%.Deblur.homd.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.homd.uncorrected.tsv Output/Metadata/%.tsv | Output/Step33 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/Macrogen.Cervix.tsv Output/Step33/Macrogen.Mouth.tsv Output/Step33/Macrogen.Neonate-3day.tsv Output/Step33/Macrogen.Neonate-5day.tsv Output/Step33/Macrogen.Vagina.tsv: Python/step33.py Output/TSV/Macrogen.tsv Output/Metadata/everything.tsv | Output/Step33 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step33/%.in: Output/Step33/%.tsv | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse format_input.py $(addprefix /,$^ $@) -c 2 -u 1 1> $@.stdout 2> $@.stderr

Output/Step33/%.res: Output/Step33/%.in | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse run_lefse.py $(addprefix /,$^ $@) --verbose 1 -y 1 1> $@.stdout 2> $@.stderr

Output/Step33/Default/%.pdf: Output/Step33/%.res | build.log Output/Step33/Default
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse plot_res.py $(addprefix /,$^ $@) --format pdf 1> $@.stdout 2> $@.stderr || $(DOCKER) touch $(addprefix /,$@)

Output/Step33/Cladogram/%.res: Output/Step33/%.res | build.log Output/Step33/Cladogram
	grep --perl-regexp "True|False" $^ > $@ || cp -f $^ $@

Output/Step33/Cladogram/%.pdf: Output/Step33/Cladogram/%.res | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse plot_cladogram.py $(addprefix /,$^ $@) --format pdf 1> $@.stdout 2> $@.stderr || $(DOCKER) touch $(addprefix /,$@)

step33-tsv: Output/Step33/everything.DADA2.gg.Cervix.tsv Output/Step33/everything.DADA2.silva.Cervix.tsv Output/Step33/everything.DADA2.homd.Cervix.tsv Output/Step33/everything.Deblur.gg.Cervix.tsv Output/Step33/everything.Deblur.silva.Cervix.tsv Output/Step33/Macrogen.Cervix.tsv
step33-in: $(addprefix Output/Step33/everything.DADA2.gg.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.DADA2.silva.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.DADA2.homd.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.Deblur.gg.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.Deblur.silva.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.Deblur.homd.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) #$(addprefix Output/Step33/Macrogen.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step33-res: $(addprefix Output/Step33/everything.DADA2.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.DADA2.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.DADA2.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.Deblur.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.Deblur.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/everything.Deblur.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) #$(addprefix Output/Step33/Macrogen.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step33-default: $(addprefix Output/Step33/Default/everything.DADA2.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Default/everything.DADA2.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Default/everything.DADA2.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Default/everything.Deblur.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Default/everything.Deblur.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Default/everything.Deblur.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) #$(addprefix Output/Step33/Default/Macrogen.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step33-cladogram: $(addprefix Output/Step33/Cladogram/everything.DADA2.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.DADA2.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.DADA2.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.Deblur.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.Deblur.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.Deblur.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.DADA2.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.DADA2.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.DADA2.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.Deblur.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.Deblur.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step33/Cladogram/everything.Deblur.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) #$(addprefix Output/Step33/Cladogram/Macrogen.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step33: step33-tsv step33-in step33-res step33-default step33-cladogram
.PHONY: step33-tsv step33-in step33-res step33-default step33-cladogram step33

# Step 34 (LEfSe analysis - Obesity)
Output/Step34 Output/Step34/Default Output/Step34/Cladogram:
	mkdir -p $@
	chmod 777 $@

Output/Step34/%.DADA2.gg.Cervix.tsv Output/Step34/%.DADA2.gg.Mouth.tsv Output/Step34/%.DADA2.gg.Neonate-3day.tsv Output/Step34/%.DADA2.gg.Neonate-5day.tsv Output/Step34/%.DADA2.gg.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.gg.uncorrected.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.DADA2.silva.Cervix.tsv Output/Step34/%.DADA2.silva.Mouth.tsv Output/Step34/%.DADA2.silva.Neonate-3day.tsv Output/Step34/%.DADA2.silva.Neonate-5day.tsv Output/Step34/%.DADA2.silva.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.silva.uncorrected.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.DADA2.homd.Cervix.tsv Output/Step34/%.DADA2.homd.Mouth.tsv Output/Step34/%.DADA2.homd.Neonate-3day.tsv Output/Step34/%.DADA2.homd.Neonate-5day.tsv Output/Step34/%.DADA2.homd.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.homd.uncorrected.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.Deblur.gg.Cervix.tsv Output/Step34/%.Deblur.gg.Mouth.tsv Output/Step34/%.Deblur.gg.Neonate-3day.tsv Output/Step34/%.Deblur.gg.Neonate-5day.tsv Output/Step34/%.Deblur.gg.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.gg.uncorrected.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.Deblur.silva.Cervix.tsv Output/Step34/%.Deblur.silva.Mouth.tsv Output/Step34/%.Deblur.silva.Neonate-3day.tsv Output/Step34/%.Deblur.silva.Neonate-5day.tsv Output/Step34/%.Deblur.silva.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.silva.uncorrected.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.Deblur.homd.Cervix.tsv Output/Step34/%.Deblur.homd.Mouth.tsv Output/Step34/%.Deblur.homd.Neonate-3day.tsv Output/Step34/%.Deblur.homd.Neonate-5day.tsv Output/Step34/%.Deblur.homd.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.homd.uncorrected.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.Macrogen.Cervix.tsv Output/Step34/%.Macrogen.Mouth.tsv Output/Step34/%.Macrogen.Neonate-3day.tsv Output/Step34/%.Macrogen.Neonate-5day.tsv Output/Step34/%.Macrogen.Vagina.tsv: Python/step33.py Output/TSV/Macrogen.tsv Output/Metadata/%.tsv | Output/Step34 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --c "Obesity" 1> $@.stdout 2> $@.stderr

Output/Step34/%.in: Output/Step34/%.tsv | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse format_input.py $(addprefix /,$^ $@) -c 2 -u 1 1> $@.stdout 2> $@.stderr

Output/Step34/%.res: Output/Step34/%.in | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse run_lefse.py $(addprefix /,$^ $@) --verbose 1 1> $@.stdout 2> $@.stderr || touch $@

Output/Step34/Default/%.pdf: Output/Step34/%.res | build.log Output/Step34/Default
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse plot_res.py $(addprefix /,$^ $@) --format pdf 1> $@.stdout 2> $@.stderr || touch $@

Output/Step34/Cladogram/%.res: Output/Step34/%.res | build.log Output/Step34/Cladogram
	grep --perl-regexp "True|False" $^ > $@ || cp -f $^ $@

Output/Step34/Cladogram/%.pdf: Output/Step34/Cladogram/%.res | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse plot_cladogram.py $(addprefix /,$^ $@) --format pdf 1> $@.stdout 2> $@.stderr || touch $@

step34-tsv: Output/Step34/everything.DADA2.gg.Cervix.tsv Output/Step34/everything.DADA2.silva.Cervix.tsv Output/Step34/everything.DADA2.homd.Cervix.tsv Output/Step34/everything.Deblur.gg.Cervix.tsv Output/Step34/everything.Deblur.silva.Cervix.tsv Output/Step34/everything.Deblur.homd.Cervix.tsv Output/Step34/everything.Macrogen.Cervix.tsv
step34-in: $(addprefix Output/Step34/everything.DADA2.gg.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.DADA2.silva.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.DADA2.homd.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Deblur.gg.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Deblur.silva.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Deblur.homd.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Macrogen.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step34-res: $(addprefix Output/Step34/everything.DADA2.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.DADA2.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.DADA2.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Deblur.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Deblur.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Deblur.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/everything.Macrogen.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step34-default: $(addprefix Output/Step34/Default/everything.DADA2.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Default/everything.DADA2.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Default/everything.DADA2.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Default/everything.Deblur.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Default/everything.Deblur.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Default/everything.Deblur.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Default/everything.Macrogen.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step34-cladogram: $(addprefix Output/Step34/Cladogram/everything.DADA2.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.DADA2.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.DADA2.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Deblur.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Deblur.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Deblur.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.DADA2.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.DADA2.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.DADA2.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Deblur.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Deblur.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Deblur.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step34/Cladogram/everything.Macrogen.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step34: step34-tsv step34-in step34-res step34-default step34-cladogram
.PHONY: step34-tsv step34-in step34-res step34-default step34-cladogram step34

# Step 35 (Clustermap plot within bacteria)
Output/Step35:
	mkdir -p $@

Output/Step35/%.pearson.png: Python/step35.py Output/TSV/%.tar.gz | build.log Output/Step35
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step35/%.spearman.png: Python/step35.py Output/TSV/%.tar.gz | build.log Output/Step35
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step35/%.kendall.png: Python/step35.py Output/TSV/%.tar.gz | build.log Output/Step35
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

step35-everything-DADA2: Output/Step35/everything.DADA2.gg.pearson.png Output/Step35/everything.DADA2.gg.spearman.png Output/Step35/everything.DADA2.gg.kendall.png Output/Step35/everything.DADA2.silva.pearson.png Output/Step35/everything.DADA2.silva.spearman.png Output/Step35/everything.DADA2.silva.kendall.png Output/Step35/everything.DADA2.homd.pearson.png Output/Step35/everything.DADA2.homd.spearman.png Output/Step35/everything.DADA2.homd.kendall.png
step35-everything-Deblur: Output/Step35/everything.Deblur.gg.pearson.png Output/Step35/everything.Deblur.gg.spearman.png Output/Step35/everything.Deblur.gg.kendall.png Output/Step35/everything.Deblur.silva.pearson.png Output/Step35/everything.Deblur.silva.spearman.png Output/Step35/everything.Deblur.silva.kendall.png Output/Step35/everything.Deblur.homd.pearson.png Output/Step35/everything.Deblur.homd.spearman.png Output/Step35/everything.Deblur.homd.kendall.png
step35-Stool-DADA2: Output/Step35/Stool.DADA2.gg.pearson.png Output/Step35/Stool.DADA2.gg.spearman.png Output/Step35/Stool.DADA2.gg.kendall.png Output/Step35/Stool.DADA2.silva.pearson.png Output/Step35/Stool.DADA2.silva.spearman.png Output/Step35/Stool.DADA2.silva.kendall.png Output/Step35/Stool.DADA2.homd.pearson.png Output/Step35/Stool.DADA2.homd.spearman.png Output/Step35/Stool.DADA2.homd.kendall.png
step35-Stool-Deblur: Output/Step35/Stool.Deblur.gg.pearson.png Output/Step35/Stool.Deblur.gg.spearman.png Output/Step35/Stool.Deblur.gg.kendall.png Output/Step35/Stool.Deblur.silva.pearson.png Output/Step35/Stool.Deblur.silva.spearman.png Output/Step35/Stool.Deblur.silva.kendall.png Output/Step35/Stool.Deblur.homd.pearson.png Output/Step35/Stool.Deblur.homd.spearman.png Output/Step35/Stool.Deblur.homd.kendall.png
step35: step35-everything-DADA2 step35-everything-Deblur
.PHONY: step35-everything-DADA2 step35-everything-Deblur step35-Stool-DADA2 step35-Stool-Deblur step35

# Step 36 (Scatter plot heatmap within bacteria)
Output/Step36:
	mkdir -p $@

Output/Step36/%.pearson.png: Python/step36.py Output/TSV/%.tar.gz | build.log Output/Step36
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step36/%.spearman.png: Python/step36.py Output/TSV/%.tar.gz | build.log Output/Step36
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step36/%.kendall.png: Python/step36.py Output/TSV/%.tar.gz | build.log Output/Step36
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

step36-everything-DADA2: Output/Step36/everything.DADA2.gg.pearson.png Output/Step36/everything.DADA2.gg.spearman.png Output/Step36/everything.DADA2.gg.kendall.png Output/Step36/everything.DADA2.silva.pearson.png Output/Step36/everything.DADA2.silva.spearman.png Output/Step36/everything.DADA2.silva.kendall.png Output/Step36/everything.DADA2.homd.pearson.png Output/Step36/everything.DADA2.homd.spearman.png Output/Step36/everything.DADA2.homd.kendall.png
step36-everything-Deblur: Output/Step36/everything.Deblur.gg.pearson.png Output/Step36/everything.Deblur.gg.spearman.png Output/Step36/everything.Deblur.gg.kendall.png Output/Step36/everything.Deblur.silva.pearson.png Output/Step36/everything.Deblur.silva.spearman.png Output/Step36/everything.Deblur.silva.kendall.png Output/Step36/everything.Deblur.homd.pearson.png Output/Step36/everything.Deblur.homd.spearman.png Output/Step36/everything.Deblur.homd.kendall.png
step36-Stool-DADA2: Output/Step36/Stool.DADA2.gg.pearson.png Output/Step36/Stool.DADA2.gg.spearman.png Output/Step36/Stool.DADA2.gg.kendall.png Output/Step36/Stool.DADA2.silva.pearson.png Output/Step36/Stool.DADA2.silva.spearman.png Output/Step36/Stool.DADA2.silva.kendall.png Output/Step36/Stool.DADA2.homd.pearson.png Output/Step36/Stool.DADA2.homd.spearman.png Output/Step36/Stool.DADA2.homd.kendall.png
step36-Stool-Deblur: Output/Step36/Stool.Deblur.gg.pearson.png Output/Step36/Stool.Deblur.gg.spearman.png Output/Step36/Stool.Deblur.gg.kendall.png Output/Step36/Stool.Deblur.silva.pearson.png Output/Step36/Stool.Deblur.silva.spearman.png Output/Step36/Stool.Deblur.silva.kendall.png Output/Step36/Stool.Deblur.homd.pearson.png Output/Step36/Stool.Deblur.homd.spearman.png Output/Step36/Stool.Deblur.homd.kendall.png
step36: step36-everything-DADA2 step36-everything-Deblur
.PHONY: step36-everything-DADA2 step36-everything-Deblur step36-Stool-DADA2 step36-Stool-Deblur step36

# Step37 (Correlation within sites)
Output/Step37:
	mkdir -p $@

Output/Step37/%.pearson.pdf: Python/step37.py Output/TSV/%.tar.gz | build.log Output/Step37
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step37/%.spearman.pdf: Python/step37.py Output/TSV/%.tar.gz | build.log Output/Step37
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step37/%.kendall.pdf: Python/step37.py Output/TSV/%.tar.gz | build.log Output/Step37
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

Output/Step37/Macrogen.%.pdf: Python/step37-1.py Output/TSV/Macrogen.tar.gz | build.log Output/Step37
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --$* --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step37-everything-DADA2: Output/Step37/everything.DADA2.gg.pearson.pdf Output/Step37/everything.DADA2.gg.spearman.pdf Output/Step37/everything.DADA2.gg.kendall.pdf Output/Step37/everything.DADA2.silva.pearson.pdf Output/Step37/everything.DADA2.silva.spearman.pdf Output/Step37/everything.DADA2.silva.kendall.pdf Output/Step37/everything.DADA2.homd.pearson.pdf Output/Step37/everything.DADA2.homd.spearman.pdf Output/Step37/everything.DADA2.homd.kendall.pdf
step37-everything-Deblur: Output/Step37/everything.Deblur.gg.pearson.pdf Output/Step37/everything.Deblur.gg.spearman.pdf Output/Step37/everything.Deblur.gg.kendall.pdf Output/Step37/everything.Deblur.silva.pearson.pdf Output/Step37/everything.Deblur.silva.spearman.pdf Output/Step37/everything.Deblur.silva.kendall.pdf Output/Step37/everything.Deblur.homd.pearson.pdf Output/Step37/everything.Deblur.homd.spearman.pdf Output/Step37/everything.Deblur.homd.kendall.pdf
step37-Stool-DADA2: Output/Step37/Stool.DADA2.gg.pearson.pdf Output/Step37/Stool.DADA2.gg.spearman.pdf Output/Step37/Stool.DADA2.gg.kendall.pdf Output/Step37/Stool.DADA2.silva.pearson.pdf Output/Step37/Stool.DADA2.silva.spearman.pdf Output/Step37/Stool.DADA2.silva.kendall.pdf Output/Step37/Stool.DADA2.homd.pearson.pdf Output/Step37/Stool.DADA2.homd.spearman.pdf Output/Step37/Stool.DADA2.homd.kendall.pdf
step37-Stool-Deblur: Output/Step37/Stool.Deblur.gg.pearson.pdf Output/Step37/Stool.Deblur.gg.spearman.pdf Output/Step37/Stool.Deblur.gg.kendall.pdf Output/Step37/Stool.Deblur.silva.pearson.pdf Output/Step37/Stool.Deblur.silva.spearman.pdf Output/Step37/Stool.Deblur.silva.kendall.pdf Output/Step37/Stool.Deblur.homd.pearson.pdf Output/Step37/Stool.Deblur.homd.spearman.pdf Output/Step37/Stool.Deblur.homd.kendall.pdf
step37: step37-everything-DADA2 step37-everything-Deblur
.PHONY: step37-everything-DADA2 step37-everything-Deblur step37-Stool-DADA2 step37-Stool-Deblur step37

# Step 38 (Clustermap plot within sites)
Output/Step38:
	mkdir -p $@

Output/Step38/everything.%.png: Python/step38.py Output/TSV/everything.%.tar.gz | build.log Output/Step38
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step38/Stool.%.png: Python/step38.py Output/TSV/Stool.%.tar.gz | build.log Output/Step38
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step38/Macrogen.png: Python/step38-1.py Output/TSV/Macrogen.tar.gz | build.log Output/Step38
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step38-everything: Output/Step38/everything.DADA2.gg.png Output/Step38/everything.DADA2.silva.png Output/Step38/everything.DADA2.homd.png Output/Step38/everything.Deblur.gg.png Output/Step38/everything.Deblur.silva.png Output/Step38/everything.Deblur.homd.png
step38-Stool: Output/Step38/Stool.DADA2.gg.png Output/Step38/Stool.DADA2.silva.png Output/Step38/Stool.DADA2.homd.png Output/Step38/Stool.Deblur.gg.png Output/Step38/Stool.Deblur.silva.png Output/Step38/Stool.Deblur.homd.png
step38: step38-everything
.PHONY: step38-everything step38-Stool step38

# Step 39 (Beta-diversity heatmap plot)
Output/Step39:
	mkdir -p $@

Output/Step39/%.DADA2.bray_curtis.pdf: Python/step39.py Output/Step32/DADA2/%.bray_curtis.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/%.DADA2.jaccard.pdf: Python/step39.py Output/Step32/DADA2/%.jaccard.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/%.DADA2.unweighted_unifrac.pdf: Python/step39.py Output/Step32/DADA2/%.unweighted_unifrac.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/%.DADA2.weighted_unifrac.pdf: Python/step39.py Output/Step32/DADA2/%.weighted_unifrac.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/%.Deblur.bray_curtis.pdf: Python/step39.py Output/Step32/Deblur/%.bray_curtis.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/%.Deblur.jaccard.pdf: Python/step39.py Output/Step32/Deblur/%.jaccard.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/%.Deblur.unweighted_unifrac.pdf: Python/step39.py Output/Step32/Deblur/%.unweighted_unifrac.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/%.Deblur.weighted_unifrac.pdf: Python/step39.py Output/Step32/Deblur/%.weighted_unifrac.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/Macrogen.%.pdf: Python/step39-1.py Output/TSV/Macrogen.tar.gz Output/Metadata/everything.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --beta "$*" 1> $@.stdout 2> $@.stderr

Output/Step39/%.DADA2.bray_curtis.Sites.tar: Python/step39-2.py Output/Step32/DADA2/%.bray_curtis.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 22> $@.stderr

Output/Step39/%.DADA2.jaccard.Sites.tar: Python/step39-2.py Output/Step32/DADA2/%.jaccard.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/%.DADA2.unweighted_unifrac.Sites.tar: Python/step39-2.py Output/Step32/DADA2/%.unweighted_unifrac.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/%.DADA2.weighted_unifrac.Sites.tar: Python/step39-2.py Output/Step32/DADA2/%.weighted_unifrac.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/%.Deblur.bray_curtis.Sites.tar: Python/step39-2.py Output/Step32/Deblur/%.bray_curtis.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 22> $@.stderr

Output/Step39/%.Deblur.jaccard.Sites.tar: Python/step39-2.py Output/Step32/Deblur/%.jaccard.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/%.Deblur.unweighted_unifrac.Sites.tar: Python/step39-2.py Output/Step32/Deblur/%.unweighted_unifrac.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/%.Deblur.weighted_unifrac.Sites.tar: Python/step39-2.py Output/Step32/Deblur/%.weighted_unifrac.tsv Output/Metadata/%.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step39/Macrogen.%.Sites.tar: Python/step39-3.py Output/TSV/Macrogen.tar.gz Output/Metadata/everything.tsv | build.log Output/Step39
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --beta "$*" 1> $@.stdout 2> $@.stderr

step39-DADA2: Output/Step39/everything.DADA2.bray_curtis.pdf Output/Step39/everything.DADA2.jaccard.pdf Output/Step39/everything.DADA2.unweighted_unifrac.pdf Output/Step39/everything.DADA2.weighted_unifrac.pdf Output/Step39/everything.DADA2.bray_curtis.Sites.tar Output/Step39/everything.DADA2.jaccard.Sites.tar Output/Step39/everything.DADA2.unweighted_unifrac.Sites.tar Output/Step39/everything.DADA2.weighted_unifrac.Sites.tar
step39-Deblur: Output/Step39/everything.Deblur.bray_curtis.pdf Output/Step39/everything.Deblur.jaccard.pdf Output/Step39/everything.Deblur.unweighted_unifrac.pdf Output/Step39/everything.Deblur.weighted_unifrac.pdf Output/Step39/everything.Deblur.bray_curtis.Sites.tar Output/Step39/everything.Deblur.jaccard.Sites.tar Output/Step39/everything.Deblur.unweighted_unifrac.Sites.tar Output/Step39/everything.Deblur.weighted_unifrac.Sites.tar
step39: step39-DADA2 step39-Deblur
.PHONY: step39-DADA2 step39-Deblur step39

# Step 40 (Violin Plot for shared taxonomy)
Output/Step40:
	mkdir -p $@

Output/Step40/everything.%.pdf: Python/step40.py Output/TSV/everything.%.tar.gz | build.log Output/Step40
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step40/Macrogen.pdf: Python/step40-1.py Output/TSV/Macrogen.tar.gz | build.log Output/Step40
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step40-everything: Output/Step40/everything.DADA2.gg.pdf Output/Step40/everything.DADA2.silva.pdf Output/Step40/everything.DADA2.homd.pdf Output/Step40/everything.Deblur.gg.pdf Output/Step40/everything.Deblur.silva.pdf Output/Step40/everything.Deblur.homd.pdf
step40: step40-everything
.PHONY: step40-everything step40

# Step 41 (Heatmap Plot for correlation between Bacteria & Clinical data)
Output/Step41:
	mkdir -p $@

Output/Step41/everything.%.pearson.pdf: Python/step41.py Output/TSV/everything.%.tar.gz Output/Metadata/everything.tsv | build.log Output/Step41
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step41/everything.%.spearman.pdf: Python/step41.py Output/TSV/everything.%.tar.gz Output/Metadata/everything.tsv | build.log Output/Step41
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step41/everything.%.kendall.pdf: Python/step41.py Output/TSV/everything.%.tar.gz Output/Metadata/everything.tsv | build.log Output/Step41
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

Output/Step41/everything.%.pearson.Sites.tar: Python/step41-2.py Output/TSV/everything.%.tar.gz Output/Metadata/everything.tsv | build.log Output/Step41
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step41/everything.%.spearman.Sites.tar: Python/step41-2.py Output/TSV/everything.%.tar.gz Output/Metadata/everything.tsv | build.log Output/Step41
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step41/everything.%.kendall.Sites.tar: Python/step41-2.py Output/TSV/everything.%.tar.gz Output/Metadata/everything.tsv | build.log Output/Step41
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

Output/Step41/Stool.%.pearson.pdf: Python/step41.py Output/TSV/Stool.%.tar.gz Output/Metadata/Stool.tsv | build.log Output/Step41
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step41/Stool.%.spearman.pdf: Python/step41.py Output/TSV/Stool.%.tar.gz Output/Metadata/Stool.tsv | build.log Output/Step41
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step41/Stool.%.kendall.pdf: Python/step41.py Output/TSV/Stool.%.tar.gz Output/Metadata/Stool.tsv | build.log Output/Step41
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

Output/Step41/Macrogen.%.pdf: Python/step41-1.py Output/TSV/Macrogen.tar.gz Output/Metadata/everything.tsv | build.log Output/Step41
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --$* 1> $@.stdout 2> $@.stderr

Output/Step41/Macrogen.%.Sites.tar: Python/step41-3.py Output/TSV/Macrogen.tar.gz Output/Metadata/everything.tsv | build.log Output/Step41
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --$* 1> $@.stdout 2> $@.stderr

step41-everything-DADA2: Output/Step41/everything.DADA2.gg.pearson.pdf Output/Step41/everything.DADA2.gg.spearman.pdf Output/Step41/everything.DADA2.gg.kendall.pdf Output/Step41/everything.DADA2.silva.pearson.pdf Output/Step41/everything.DADA2.silva.spearman.pdf Output/Step41/everything.DADA2.silva.kendall.pdf Output/Step41/everything.DADA2.homd.pearson.pdf Output/Step41/everything.DADA2.homd.spearman.pdf Output/Step41/everything.DADA2.homd.kendall.pdf Output/Step41/everything.DADA2.gg.pearson.Sites.tar Output/Step41/everything.DADA2.gg.spearman.Sites.tar Output/Step41/everything.DADA2.gg.kendall.Sites.tar Output/Step41/everything.DADA2.silva.pearson.Sites.tar Output/Step41/everything.DADA2.silva.spearman.Sites.tar Output/Step41/everything.DADA2.silva.kendall.Sites.tar Output/Step41/everything.DADA2.homd.pearson.Sites.tar Output/Step41/everything.DADA2.homd.spearman.Sites.tar Output/Step41/everything.DADA2.homd.kendall.Sites.tar
step41-everything-Deblur: Output/Step41/everything.Deblur.gg.pearson.pdf Output/Step41/everything.Deblur.gg.spearman.pdf Output/Step41/everything.Deblur.gg.kendall.pdf Output/Step41/everything.Deblur.silva.pearson.pdf Output/Step41/everything.Deblur.silva.spearman.pdf Output/Step41/everything.Deblur.silva.kendall.pdf Output/Step41/everything.Deblur.homd.pearson.pdf Output/Step41/everything.Deblur.homd.spearman.pdf Output/Step41/everything.Deblur.homd.kendall.pdf Output/Step41/everything.Deblur.gg.pearson.Sites.tar Output/Step41/everything.Deblur.gg.spearman.Sites.tar Output/Step41/everything.Deblur.gg.kendall.Sites.tar Output/Step41/everything.Deblur.silva.pearson.Sites.tar Output/Step41/everything.Deblur.silva.spearman.Sites.tar Output/Step41/everything.Deblur.silva.kendall.Sites.tar Output/Step41/everything.Deblur.homd.pearson.Sites.tar Output/Step41/everything.Deblur.homd.spearman.Sites.tar Output/Step41/everything.Deblur.homd.kendall.Sites.tar
step41: step41-everything-DADA2 step41-everything-Deblur
.PHONY: step41-everything-DADA2 step41-everything-Deblur step41

# Step 42 (Violin plot for taxonomy with clinical data)
Output/Step42:
	mkdir -p $@

Output/Step42/everything.%.tar: Python/step42.py Output/TSV/everything.%.tar.gz Output/Metadata/everything.tsv | Output/Step42 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step42-DADA2: Output/Step42/everything.DADA2.gg.tar Output/Step42/everything.DADA2.silva.tar Output/Step42/everything.DADA2.homd.tar
step42-Deblur: Output/Step42/everything.Deblur.gg.tar Output/Step42/everything.Deblur.silva.tar Output/Step42/everything.Deblur.homd.tar
step42: step42-DADA2 step42-Deblur
.PHONY: step42-DADA2 step42-Deblur step42

# Step 43 (Volcano plot with Corrleation)
Output/Step43:
	mkdir -p $@

Output/Step43/everything.%.pearson.tar: Python/step43.py Output/TSV/everything.%.tar.gz Output/Metadata/everything.tsv | build.log Output/Step43
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step43/everything.%.spearman.tar: Python/step43.py Output/TSV/everything.%.tar.gz Output/Metadata/everything.tsv | build.log Output/Step43
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step43/everything.%.kendall.tar: Python/step43.py Output/TSV/everything.%.tar.gz Output/Metadata/everything.tsv | build.log Output/Step43
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

step43-everything-DADA2: Output/Step43/everything.DADA2.gg.pearson.tar Output/Step43/everything.DADA2.silva.pearson.tar Output/Step43/everything.DADA2.homd.pearson.tar Output/Step43/everything.DADA2.gg.spearman.tar Output/Step43/everything.DADA2.silva.spearman.tar Output/Step43/everything.DADA2.homd.spearman.tar Output/Step43/everything.DADA2.gg.kendall.tar Output/Step43/everything.DADA2.silva.kendall.tar Output/Step43/everything.DADA2.homd.kendall.tar
step43-everything-Deblur: Output/Step43/everything.Deblur.gg.pearson.tar Output/Step43/everything.Deblur.silva.pearson.tar Output/Step43/everything.Deblur.homd.pearson.tar Output/Step43/everything.Deblur.gg.spearman.tar Output/Step43/everything.Deblur.silva.spearman.tar Output/Step43/everything.Deblur.homd.spearman.tar Output/Step43/everything.Deblur.gg.kendall.tar Output/Step43/everything.Deblur.silva.kendall.tar Output/Step43/everything.Deblur.homd.kendall.tar
step43: step43-everything-DADA2 step43-everything-Deblur
.PHONY: step43-everything-DADA2 step43-everything-Deblur step43

# Step 44 (Volcano plot with differetially abundant taxa)
Output/Step44 Output/Step44_Proportion:
	mkdir -p $@

Output/Step44/everything.%.Cervix.tsv Output/Step44/everything.%.Mouth.tsv Output/Step44/everything.%.Vagina.tsv Output/Step44/everything.%.Cervix.coldata Output/Step44/everything.%.Mouth.coldata Output/Step44/everything.%.Vagina.coldata: Python/step44-2.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step44
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step44/everything.%.EL.Cervix.tsv Output/Step44/everything.%.EL.Mouth.tsv Output/Step44/everything.%.EL.Vagina.tsv Output/Step44/everything.%.EL.Cervix.coldata Output/Step44/everything.%.EL.Mouth.coldata Output/Step44/everything.%.EL.Vagina.coldata: Python/step44-2.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step44
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --compare "Detail Premature" "Late PTB" "Early PTB" 1> $@.stdout 2> $@.stderr

Output/Step44/everything.%.EF.Cervix.tsv Output/Step44/everything.%.EF.Mouth.tsv Output/Step44/everything.%.EF.Vagina.tsv Output/Step44/everything.%.EF.Cervix.coldata Output/Step44/everything.%.EF.Mouth.coldata Output/Step44/everything.%.EF.Vagina.coldata: Python/step44-2.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step44
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --compare "Detail Premature" "Normal" "Early PTB" 1> $@.stdout 2> $@.stderr

Output/Step44/everything.%.LF.Cervix.tsv Output/Step44/everything.%.LF.Mouth.tsv Output/Step44/everything.%.LF.Vagina.tsv Output/Step44/everything.%.LF.Cervix.coldata Output/Step44/everything.%.LF.Mouth.coldata Output/Step44/everything.%.LF.Vagina.coldata: Python/step44-2.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step44
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --compare "Detail Premature" "Normal" "Late PTB" 1> $@.stdout 2> $@.stderr

Output/Step44/everything.%.Cervix.DAT.tsv: R/run_DESeq.R Output/Step44/everything.%.Cervix.tsv Output/Step44/everything.%.Cervix.coldata | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --readcount $(addprefix /,$(filter %.tsv,$^)) --coldata $(addprefix /,$(filter %.coldata,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step44/everything.%.Mouth.DAT.tsv: R/run_DESeq.R Output/Step44/everything.%.Mouth.tsv Output/Step44/everything.%.Mouth.coldata | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --readcount $(addprefix /,$(filter %.tsv,$^)) --coldata $(addprefix /,$(filter %.coldata,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step44/everything.%.Vagina.DAT.tsv: R/run_DESeq.R Output/Step44/everything.%.Vagina.tsv Output/Step44/everything.%.Vagina.coldata | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --readcount $(addprefix /,$(filter %.tsv,$^)) --coldata $(addprefix /,$(filter %.coldata,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step44/everything.%.pdf: Python/step44.py Output/Step44/everything.%.DAT.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step44_Proportion/everything.%.Cervix.tsv Output/Step44_Proportion/everything.%.Mouth.tsv Output/Step44_Proportion/everything.%.Vagina.tsv Output/Step44_Proportion/everything.%.Cervix.coldata Output/Step44_Proportion/everything.%.Mouth.coldata Output/Step44_Proportion/everything.%.Vagina.coldata: Python/step44-2.py Output/TSV_Proportion/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step44_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step44_Proportion/everything.%.EL.Cervix.tsv Output/Step44_Proportion/everything.%.EL.Mouth.tsv Output/Step44_Proportion/everything.%.EL.Vagina.tsv Output/Step44_Proportion/everything.%.EL.Cervix.coldata Output/Step44_Proportion/everything.%.EL.Mouth.coldata Output/Step44_Proportion/everything.%.EL.Vagina.coldata: Python/step44-2.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step44_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --compare "Detail Premature" "Late PTB" "Early PTB" 1> $@.stdout 2> $@.stderr

Output/Step44_Proportion/everything.%.EF.Cervix.tsv Output/Step44_Proportion/everything.%.EF.Mouth.tsv Output/Step44_Proportion/everything.%.EF.Vagina.tsv Output/Step44_Proportion/everything.%.EF.Cervix.coldata Output/Step44_Proportion/everything.%.EF.Mouth.coldata Output/Step44_Proportion/everything.%.EF.Vagina.coldata: Python/step44-2.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step44_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --compare "Detail Premature" "Normal" "Early PTB" 1> $@.stdout 2> $@.stderr

Output/Step44_Proportion/everything.%.LF.Cervix.tsv Output/Step44_Proportion/everything.%.LF.Mouth.tsv Output/Step44_Proportion/everything.%.LF.Vagina.tsv Output/Step44_Proportion/everything.%.LF.Cervix.coldata Output/Step44_Proportion/everything.%.LF.Mouth.coldata Output/Step44_Proportion/everything.%.LF.Vagina.coldata: Python/step44-2.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step44_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) --compare "Detail Premature" "Normal" "Late PTB" 1> $@.stdout 2> $@.stderr

Output/Step44_Proportion/everything.%.Cervix.DAT.tsv: R/run_DESeq.R Output/Step44_Proportion/everything.%.Cervix.tsv Output/Step44_Proportion/everything.%.Cervix.coldata | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --readcount $(addprefix /,$(filter %.tsv,$^)) --coldata $(addprefix /,$(filter %.coldata,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step44_Proportion/everything.%.Mouth.DAT.tsv: R/run_DESeq.R Output/Step44_Proportion/everything.%.Mouth.tsv Output/Step44_Proportion/everything.%.Mouth.coldata | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --readcount $(addprefix /,$(filter %.tsv,$^)) --coldata $(addprefix /,$(filter %.coldata,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step44_Proportion/everything.%.Vagina.DAT.tsv: R/run_DESeq.R Output/Step44_Proportion/everything.%.Vagina.tsv Output/Step44_Proportion/everything.%.Vagina.coldata | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --readcount $(addprefix /,$(filter %.tsv,$^)) --coldata $(addprefix /,$(filter %.coldata,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step44_Proportion/everything.%.pdf: Python/step44.py Output/Step44_Proportion/everything.%.DAT.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step44-uncorrected: Output/Step44/everything.DADA2.homd.uncorrected.Cervix.coldata $(addprefix Output/Step44/everything.DADA2.homd.uncorrected.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44/everything.DADA2.homd.uncorrected.,$(addsuffix .pdf,Cervix Mouth Vagina)) Output/Step44/everything.DADA2.homd.uncorrected.EL.Cervix.coldata $(addprefix Output/Step44/everything.DADA2.homd.uncorrected.EL.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44/everything.DADA2.homd.uncorrected.EL.,$(addsuffix .pdf,Cervix Mouth Vagina)) Output/Step44/everything.DADA2.homd.uncorrected.EF.Cervix.coldata $(addprefix Output/Step44/everything.DADA2.homd.uncorrected.EF.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44/everything.DADA2.homd.uncorrected.EF.,$(addsuffix .pdf,Cervix Mouth Vagina)) Output/Step44/everything.DADA2.homd.uncorrected.LF.Cervix.coldata $(addprefix Output/Step44/everything.DADA2.homd.uncorrected.LF.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44/everything.DADA2.homd.uncorrected.LF.,$(addsuffix .pdf,Cervix Mouth Vagina))
step44-corrected: Output/Step44/everything.DADA2.homd.Cervix.coldata $(addprefix Output/Step44/everything.DADA2.homd.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44/everything.DADA2.homd.,$(addsuffix .pdf,Cervix Mouth Vagina)) Output/Step44/everything.DADA2.homd.EL.Cervix.coldata $(addprefix Output/Step44/everything.DADA2.homd.EL.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44/everything.DADA2.homd.EL.,$(addsuffix .pdf,Cervix Mouth Vagina)) Output/Step44/everything.DADA2.homd.EF.Cervix.coldata $(addprefix Output/Step44/everything.DADA2.homd.EF.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44/everything.DADA2.homd.EF.,$(addsuffix .pdf,Cervix Mouth Vagina)) Output/Step44/everything.DADA2.homd.LF.Cervix.coldata $(addprefix Output/Step44/everything.DADA2.homd.LF.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44/everything.DADA2.homd.LF.,$(addsuffix .pdf,Cervix Mouth Vagina))
step44-uncorrected-Proportion: Output/Step44_Proportion/everything.DADA2.homd.uncorrected.Cervix.coldata $(addprefix Output/Step44_Proportion/everything.DADA2.homd.uncorrected.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44_Proportion/everything.DADA2.homd.uncorrected.,$(addsuffix .pdf,Cervix Mouth Vagina)) Output/Step44_Proportion/everything.DADA2.homd.uncorrected.EL.Cervix.coldata $(addprefix Output/Step44_Proportion/everything.DADA2.homd.uncorrected.EL.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44_Proportion/everything.DADA2.homd.uncorrected.EL.,$(addsuffix .pdf,Cervix Mouth Vagina)) Output/Step44_Proportion/everything.DADA2.homd.uncorrected.EF.Cervix.coldata $(addprefix Output/Step44_Proportion/everything.DADA2.homd.uncorrected.EF.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44_Proportion/everything.DADA2.homd.uncorrected.EF.,$(addsuffix .pdf,Cervix Mouth Vagina)) Output/Step44_Proportion/everything.DADA2.homd.uncorrected.LF.Cervix.coldata $(addprefix Output/Step44_Proportion/everything.DADA2.homd.uncorrected.LF.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44_Proportion/everything.DADA2.homd.uncorrected.LF.,$(addsuffix .pdf,Cervix Mouth Vagina))
step44-corrected-Proportion: Output/Step44_Proportion/everything.DADA2.homd.Cervix.coldata $(addprefix Output/Step44_Proportion/everything.DADA2.homd.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44_Proportion/everything.DADA2.homd.,$(addsuffix .pdf,Cervix Mouth Vagina)) Output/Step44_Proportion/everything.DADA2.homd.EL.Cervix.coldata $(addprefix Output/Step44_Proportion/everything.DADA2.homd.EL.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44_Proportion/everything.DADA2.homd.EL.,$(addsuffix .pdf,Cervix Mouth Vagina)) Output/Step44_Proportion/everything.DADA2.homd.EF.Cervix.coldata $(addprefix Output/Step44_Proportion/everything.DADA2.homd.EF.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44_Proportion/everything.DADA2.homd.EF.,$(addsuffix .pdf,Cervix Mouth Vagina)) Output/Step44_Proportion/everything.DADA2.homd.LF.Cervix.coldata $(addprefix Output/Step44_Proportion/everything.DADA2.homd.LF.,$(addsuffix .DAT.tsv,Cervix Mouth Vagina)) $(addprefix Output/Step44_Proportion/everything.DADA2.homd.LF.,$(addsuffix .pdf,Cervix Mouth Vagina))
step44: step44-uncorrected step44-corrected
.PHONY: step44-uncorrected step44-corrected step44

# Step 45 (Alpha-diversity)
Output/Step45:
	mkdir -p $@

Output/Step45/%/tree.nwk: Output/Step28/%.rooted-tree.qza | Output/Step45 build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$^) --output-path $(addprefix /,$(@D)) 1> $(@D).stdout 2> $(@D).stderr

Output/Step45/%.uncorrected.tar: Python/step45.py Output/TSV/%.homd.uncorrected.tsv Output/Metadata/everything.tsv Output/Step45/%/tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step45/%.tar: Python/step45.py Output/TSV/%.homd.tsv Output/Metadata/everything.tsv Output/Step45/%/tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step45/%.uncorrected.tsv: Python/step45-1.py Output/TSV/%.homd.uncorrected.tsv Output/Metadata/everything.tsv Output/Step45/%/tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step45/%.tsv: Python/step45-1.py Output/TSV/%.homd.tsv Output/Metadata/everything.tsv Output/Step45/%/tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step45-tree: Output/Step45/everything.DADA2/tree.nwk Output/Step45/everything.Deblur/tree.nwk
step45-alpha: Output/Step45/everything.DADA2.uncorrected.tar Output/Step45/everything.DADA2.tar Output/Step45/everything.Deblur.uncorrected.tar Output/Step45/everything.Deblur.tar
step45-alpha-tsv: Output/Step45/everything.DADA2.uncorrected.tsv Output/Step45/everything.DADA2.tsv Output/Step45/everything.Deblur.uncorrected.tsv Output/Step45/everything.Deblur.tsv
step45: step45-tree step45-alpha step45-alpha-tsv
.PHONY: step45-tree step45-alpha step45-alpha-tsv step45

# Step 46 (Beta-diversity)
Output/Step46:
	mkdir -p $@

Output/Step46/%/tree.nwk: Output/Step28/%.rooted-tree.qza | Output/Step46 build.log
	$(DOCKER) qiime tools export --input-path $(addprefix /,$^) --output-path $(addprefix /,$(@D)) 1> $(@D).stdout 2> $(@D).stderr

Output/Step46/%.uncorrected.tar: Python/step46.py Output/TSV/%.homd.uncorrected.tsv Output/Metadata/everything.tsv Output/Step46/%/tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step46/%.tar: Python/step46.py Output/TSV/%.homd.tsv Output/Metadata/everything.tsv Output/Step46/%/tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step46/%.uncorrected.unweighted_unifrac.tsv Output/Step46/%.uncorrected.weighted_unifrac.tsv: Python/step46-1.py Output/TSV/%.homd.uncorrected.tsv Output/Metadata/everything.tsv Output/Step46/%/tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step46/%.unweighted_unifrac.tsv Output/Step46/%.weighted_unifrac.tsv: Python/step46-1.py Output/TSV/%.homd.tsv Output/Metadata/everything.tsv Output/Step46/%/tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

step46-tree: Output/Step46/everything.DADA2/tree.nwk Output/Step46/everything.Deblur/tree.nwk
step46-beta: Output/Step46/everything.DADA2.uncorrected.tar Output/Step46/everything.DADA2.tar Output/Step46/everything.Deblur.uncorrected.tar Output/Step46/everything.Deblur.tar
step46-beta-tsv: Output/Step46/everything.DADA2.uncorrected.unweighted_unifrac.tsv Output/Step46/everything.DADA2.unweighted_unifrac.tsv Output/Step46/everything.Deblur.uncorrected.unweighted_unifrac.tsv Output/Step46/everything.Deblur.unweighted_unifrac.tsv
step46: step46-tree step46-beta step46-beta-tsv
.PHONY: step46-tree step46-beta step46-beta-tsv step46

# Step 47 (Beta-diversity heatmap plot)
Output/Step47:
	mkdir -p $@

Output/Step47/%.pdf: Python/step47.py Output/Step46/%.tsv Output/Metadata/everything.tsv | Output/Step47 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step47/%.Sites.tar: Python/step47-2.py Output/Step46/%.tsv Output/Metadata/everything.tsv | Output/Step47 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step47-All-uncorrected: Output/Step47/everything.DADA2.uncorrected.unweighted_unifrac.pdf Output/Step47/everything.DADA2.uncorrected.weighted_unifrac.pdf Output/Step47/everything.Deblur.uncorrected.unweighted_unifrac.pdf Output/Step47/everything.Deblur.uncorrected.weighted_unifrac.pdf
step47-All: Output/Step47/everything.DADA2.unweighted_unifrac.pdf Output/Step47/everything.DADA2.weighted_unifrac.pdf Output/Step47/everything.Deblur.unweighted_unifrac.pdf Output/Step47/everything.Deblur.weighted_unifrac.pdf
step47-Sites: Output/Step47/everything.DADA2.uncorrected.unweighted_unifrac.Sites.tar Output/Step47/everything.DADA2.uncorrected.weighted_unifrac.Sites.tar Output/Step47/everything.Deblur.uncorrected.unweighted_unifrac.Sites.tar Output/Step47/everything.Deblur.uncorrected.weighted_unifrac.Sites.tar
step47-Sites: Output/Step47/everything.DADA2.unweighted_unifrac.Sites.tar Output/Step47/everything.DADA2.weighted_unifrac.Sites.tar Output/Step47/everything.Deblur.unweighted_unifrac.Sites.tar Output/Step47/everything.Deblur.weighted_unifrac.Sites.tar
step47: step47-All-uncorrected step47-All step47-Sites-uncorrected step47-Sites
.PHONY: step47-All-uncorrected step47-All step47-Sites-uncorrected step47-Sites step47

# Step 48 (Violin plot for Shared taxonomy proportion)
Output/Step48:
	mkdir -p $@

Output/Step48/everything.%.pdf: Python/step48.py Output/TSV/everything.%.tsv | build.log Output/Step48
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step48: Output/Step48/everything.DADA2.gg.pdf Output/Step48/everything.DADA2.silva.pdf Output/Step48/everything.DADA2.homd.pdf Output/Step48/everything.Deblur.gg.pdf Output/Step48/everything.Deblur.silva.pdf Output/Step48/everything.Deblur.homd.pdf
.PHONY: step48

# Step 49 (Heatmap Plot for correlation between Bacteria & Clinical data)
Output/Step49 Output/Step49_Proportion:
	mkdir -p $@

Output/Step49/everything.%.pearson.pdf: Python/step49.py Output/TSV/everything.%.uncorrected.tsv Output/Metadata/everything.tsv | build.log Output/Step49
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step49/everything.%.spearman.pdf: Python/step49.py Output/TSV/everything.%.uncorrected.tsv Output/Metadata/everything.tsv | build.log Output/Step49
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step49/everything.%.kendall.pdf: Python/step49.py Output/TSV/everything.%.uncorrected.tsv Output/Metadata/everything.tsv | build.log Output/Step49
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

Output/Step49/everything.%.pearson.Sites.tar: Python/step49-2.py Output/TSV/everything.%.uncorrected.tsv Output/Metadata/everything.tsv | build.log Output/Step49
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step49/everything.%.spearman.Sites.tar: Python/step49-2.py Output/TSV/everything.%.uncorrected.tsv Output/Metadata/everything.tsv | build.log Output/Step49
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step49/everything.%.kendall.Sites.tar: Python/step49-2.py Output/TSV/everything.%.uncorrected.tsv Output/Metadata/everything.tsv | build.log Output/Step49
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

Output/Step49_Proportion/everything.%.pearson.pdf: Python/step49.py Output/TSV/everything.%.uncorrected.tsv Output/Metadata/everything.tsv | build.log Output/Step49_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step49_Proportion/everything.%.spearman.pdf: Python/step49.py Output/TSV/everything.%.uncorrected.tsv Output/Metadata/everything.tsv | build.log Output/Step49_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step49_Proportion/everything.%.kendall.pdf: Python/step49.py Output/TSV/everything.%.uncorrected.tsv Output/Metadata/everything.tsv | build.log Output/Step49_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

Output/Step49_Proportion/everything.%.pearson.Sites.tar: Python/step49-2.py Output/TSV/everything.%.uncorrected.tsv Output/Metadata/everything.tsv | build.log Output/Step49_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --pearson 1> $@.stdout 2> $@.stderr

Output/Step49_Proportion/everything.%.spearman.Sites.tar: Python/step49-2.py Output/TSV/everything.%.uncorrected.tsv Output/Metadata/everything.tsv | build.log Output/Step49_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --spearman 1> $@.stdout 2> $@.stderr

Output/Step49_Proportion/everything.%.kendall.Sites.tar: Python/step49-2.py Output/TSV/everything.%.uncorrected.tsv Output/Metadata/everything.tsv | build.log Output/Step49_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --kendall 1> $@.stdout 2> $@.stderr

step49-DADA2: Output/Step49/everything.DADA2.gg.pearson.pdf Output/Step49/everything.DADA2.gg.spearman.pdf Output/Step49/everything.DADA2.gg.kendall.pdf Output/Step49/everything.DADA2.silva.pearson.pdf Output/Step49/everything.DADA2.silva.spearman.pdf Output/Step49/everything.DADA2.silva.kendall.pdf Output/Step49/everything.DADA2.homd.pearson.pdf Output/Step49/everything.DADA2.homd.spearman.pdf Output/Step49/everything.DADA2.homd.kendall.pdf Output/Step49/everything.DADA2.gg.pearson.Sites.tar Output/Step49/everything.DADA2.gg.spearman.Sites.tar Output/Step49/everything.DADA2.gg.kendall.Sites.tar Output/Step49/everything.DADA2.silva.pearson.Sites.tar Output/Step49/everything.DADA2.silva.spearman.Sites.tar Output/Step49/everything.DADA2.silva.kendall.Sites.tar Output/Step49/everything.DADA2.homd.pearson.Sites.tar Output/Step49/everything.DADA2.homd.spearman.Sites.tar Output/Step49/everything.DADA2.homd.kendall.Sites.tar
step49-Deblur: Output/Step49/everything.Deblur.gg.pearson.pdf Output/Step49/everything.Deblur.gg.spearman.pdf Output/Step49/everything.Deblur.gg.kendall.pdf Output/Step49/everything.Deblur.silva.pearson.pdf Output/Step49/everything.Deblur.silva.spearman.pdf Output/Step49/everything.Deblur.silva.kendall.pdf Output/Step49/everything.Deblur.homd.pearson.pdf Output/Step49/everything.Deblur.homd.spearman.pdf Output/Step49/everything.Deblur.homd.kendall.pdf Output/Step49/everything.Deblur.gg.pearson.Sites.tar Output/Step49/everything.Deblur.gg.spearman.Sites.tar Output/Step49/everything.Deblur.gg.kendall.Sites.tar Output/Step49/everything.Deblur.silva.pearson.Sites.tar Output/Step49/everything.Deblur.silva.spearman.Sites.tar Output/Step49/everything.Deblur.silva.kendall.Sites.tar Output/Step49/everything.Deblur.homd.pearson.Sites.tar Output/Step49/everything.Deblur.homd.spearman.Sites.tar Output/Step49/everything.Deblur.homd.kendall.Sites.tar
step49-DADA2-Proportion: Output/Step49_Proportion/everything.DADA2.gg.pearson.pdf Output/Step49_Proportion/everything.DADA2.gg.spearman.pdf Output/Step49_Proportion/everything.DADA2.gg.kendall.pdf Output/Step49_Proportion/everything.DADA2.silva.pearson.pdf Output/Step49_Proportion/everything.DADA2.silva.spearman.pdf Output/Step49_Proportion/everything.DADA2.silva.kendall.pdf Output/Step49_Proportion/everything.DADA2.homd.pearson.pdf Output/Step49_Proportion/everything.DADA2.homd.spearman.pdf Output/Step49_Proportion/everything.DADA2.homd.kendall.pdf Output/Step49_Proportion/everything.DADA2.gg.pearson.Sites.tar Output/Step49_Proportion/everything.DADA2.gg.spearman.Sites.tar Output/Step49_Proportion/everything.DADA2.gg.kendall.Sites.tar Output/Step49_Proportion/everything.DADA2.silva.pearson.Sites.tar Output/Step49_Proportion/everything.DADA2.silva.spearman.Sites.tar Output/Step49_Proportion/everything.DADA2.silva.kendall.Sites.tar Output/Step49_Proportion/everything.DADA2.homd.pearson.Sites.tar Output/Step49_Proportion/everything.DADA2.homd.spearman.Sites.tar Output/Step49_Proportion/everything.DADA2.homd.kendall.Sites.tar
step49-Deblur-Proportion: Output/Step49_Proportion/everything.Deblur.gg.pearson.pdf Output/Step49_Proportion/everything.Deblur.gg.spearman.pdf Output/Step49_Proportion/everything.Deblur.gg.kendall.pdf Output/Step49_Proportion/everything.Deblur.silva.pearson.pdf Output/Step49_Proportion/everything.Deblur.silva.spearman.pdf Output/Step49_Proportion/everything.Deblur.silva.kendall.pdf Output/Step49_Proportion/everything.Deblur.homd.pearson.pdf Output/Step49_Proportion/everything.Deblur.homd.spearman.pdf Output/Step49_Proportion/everything.Deblur.homd.kendall.pdf Output/Step49_Proportion/everything.Deblur.gg.pearson.Sites.tar Output/Step49_Proportion/everything.Deblur.gg.spearman.Sites.tar Output/Step49_Proportion/everything.Deblur.gg.kendall.Sites.tar Output/Step49_Proportion/everything.Deblur.silva.pearson.Sites.tar Output/Step49_Proportion/everything.Deblur.silva.spearman.Sites.tar Output/Step49_Proportion/everything.Deblur.silva.kendall.Sites.tar Output/Step49_Proportion/everything.Deblur.homd.pearson.Sites.tar Output/Step49_Proportion/everything.Deblur.homd.spearman.Sites.tar Output/Step49_Proportion/everything.Deblur.homd.kendall.Sites.tar
step49: step49-DADA2 step49-Deblur step49-DADA2-Proportion step49-Deblur-Proportion
.PHONY: step49-DADA2 step49-Deblur step49-DADA2-Proportion step49-Deblur-Proportion step49

# Step 50 (Volcano plot with differentially expressed taxa)
Output/Step50:
	mkdir -p $@

Output/Step50/everything.%.Mouth.tsv Output/Step50/everything.%.Mouth.coldata Output/Step50/everything.%.Cervix.tsv Output/Step50/everything.%.Cervix.coldata Output/Step50/everything.%.Vagina.tsv Output/Step50/everything.%.Vagina.coldata Output/Step50/everything.%.Neonate-1day.tsv Output/Step50/everything.%.Neonate-1day.coldata Output/Step50/everything.%.Neonate-3day.tsv Output/Step50/everything.%.Neonate-3day.coldata Output/Step50/everything.%.Neonate-5day.tsv Output/Step50/everything.%.Neonate-5day.coldata: Python/step50-1.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step50
	$(DOCKER) $(PYTHON) $(addprefix /,$^ Output/Step50/everything.$*) 1> $@.stdout 2> $@.stdout

Output/Step50/everything.%.DEG.tsv: R/run_DESeq.R Output/Step50/everything.%.tsv Output/Step50/everything.%.coldata | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --readcount $(addprefix /,$(filter %.tsv,$^)) --coldata $(addprefix /,$(filter %.coldata,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step50/everything.%.pdf: Python/step50.py Output/Step50/everything.%.tsv Output/Metadata/everything.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step50-uncorrected: Output/Step50/everything.DADA2.homd.uncorrected.Mouth.tsv Output/Step50/everything.DADA2.homd.uncorrected.Mouth.DEG.tsv Output/Step50/everything.DADA2.homd.uncorrected.Mouth.DEG.pdf Output/Step50/everything.DADA2.homd.uncorrected.Neonate-3day.tsv Output/Step50/everything.DADA2.homd.uncorrected.Neonate-3day.DEG.tsv Output/Step50/everything.DADA2.homd.uncorrected.Neonate-3day.DEG.pdf Output/Step50/everything.DADA2.homd.uncorrected.Neonate-5day.tsv Output/Step50/everything.DADA2.homd.uncorrected.Neonate-5day.DEG.tsv Output/Step50/everything.DADA2.homd.uncorrected.Neonate-5day.DEG.pdf
step50-corrected: Output/Step50/everything.DADA2.homd.Mouth.tsv Output/Step50/everything.DADA2.homd.Mouth.DEG.tsv Output/Step50/everything.DADA2.homd.Mouth.DEG.pdf Output/Step50/everything.DADA2.homd.Neonate-3day.tsv Output/Step50/everything.DADA2.homd.Neonate-3day.DEG.tsv Output/Step50/everything.DADA2.homd.Neonate-3day.DEG.pdf Output/Step50/everything.DADA2.homd.Neonate-5day.tsv Output/Step50/everything.DADA2.homd.Neonate-5day.DEG.tsv Output/Step50/everything.DADA2.homd.Neonate-5day.DEG.pdf
step50: step50-uncorrected step50-corrected
.PHONY: step50

# Step 51 (QIIME2 - Batch corrected)
Output/Step51:
	mkdir -p $@

Output/Step51/%.biom: Output/TSV/%.tsv | build.log Output/Step51
	$(DOCKER) biom convert --input-fp $(addprefix /,$^) --output-fp $(addprefix /,$@) --table-type="OTU table" --to-hdf5

Output/Step51/%.frequency.qza: Output/Step51/%.biom | build.log
	$(DOCKER) qiime tools import --input-path $(addprefix /,$^) --type "FeatureTable[Frequency]" --input-format "BIOMV210Format" --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step51/%.tsv: Python/step51.py Output/TSV/%.tsv | build.log Output/Step51
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step51/%.taxonomy.qza: Output/Step51/%.tsv | build.log
	$(DOCKER) qiime tools import --input-path $(addprefix /,$^) --type "FeatureData[Taxonomy]" --input-format "TSVTaxonomyFormat" --output-path $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step51-biom: Output/Step51/everything.DADA2.gg.biom Output/Step51/everything.DADA2.silva.biom Output/Step51/everything.DADA2.homd.biom Output/Step51/everything.Deblur.gg.biom Output/Step51/everything.Deblur.silva.biom Output/Step51/everything.Deblur.homd.biom Output/Step51/Macrogen.biom
step51-frequency-qza: Output/Step51/everything.DADA2.gg.frequency.qza Output/Step51/everything.DADA2.silva.frequency.qza Output/Step51/everything.DADA2.homd.frequency.qza Output/Step51/everything.Deblur.gg.frequency.qza Output/Step51/everything.Deblur.silva.frequency.qza Output/Step51/everything.Deblur.homd.frequency.qza Output/Step51/Macrogen.frequency.qza
step51-taxonomy-qza: Output/Step51/everything.DADA2.gg.tsv Output/Step51/everything.DADA2.silva.tsv Output/Step51/everything.DADA2.homd.tsv Output/Step51/everything.Deblur.gg.tsv Output/Step51/everything.Deblur.silva.tsv Output/Step51/everything.Deblur.homd.tsv Output/Step51/Macrogen.tsv
step51: step51-biom step51-frequency-qza step51-taxonomy-qza
.PHONY: step51-biom step51-frequency-qza step51-taxonomy-qza step51

# Step 52 (Run ANCOM - Batch corrected)
Output/Step52:
	mkdir -p $@

Output/Step52/%.DADA2.gg.qza: Output/DADA2/%.table.qza Output/Step51/%.DADA2.gg.taxonomy.qza | build.log Output/Step52
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/%.DADA2.silva.qza: Output/DADA2/%.table.qza Output/Step51/%.DADA2.silva.taxonomy.qza | build.log Output/Step52
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/%.DADA2.homd.qza: Output/DADA2/%.table.qza Output/Step51/%.DADA2.homd.taxonomy.qza | build.log Output/Step52
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/%.Deblur.gg.qza: Output/Deblur/%.table.qza Output/Step51/%.Deblur.gg.taxonomy.qza | build.log Output/Step52
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/%.Deblur.silva.qza: Output/Deblur/%.table.qza Output/Step51/%.Deblur.silva.taxonomy.qza | build.log Output/Step52
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/%.Deblur.homd.qza: Output/Deblur/%.table.qza Output/Step51/%.Deblur.homd.taxonomy.qza | build.log Output/Step52
	$(DOCKER) qiime taxa collapse --i-table $(addprefix /,$(word 1,$^)) --i-taxonomy $(addprefix /,$(word 2,$^)) --p-level 7 --o-collapsed-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/everything.%.Mouth.qza: Output/Step52/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]='Mouth'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/everything.%.Cervix.qza: Output/Step52/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]='Cervix'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/everything.%.Vagina.qza: Output/Step52/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]='Vagina'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/everything.%.Neonate.qza: Output/Step52/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site] IN ('Neonate-1day', 'Neonate-3day', 'Neonate-5day')" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/everything.%.Neonate-1day.qza: Output/Step52/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]=='Neonate-1day'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/everything.%.Neonate-3day.qza: Output/Step52/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]=='Neonate-3day'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/everything.%.Neonate-5day.qza: Output/Step52/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime feature-table filter-samples --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --p-where "[Site]=='Neonate-5day'" --o-filtered-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/%.pseudocount.qza: Output/Step52/%.qza | build.log
	$(DOCKER) qiime composition add-pseudocount --i-table $(addprefix /,$(word 1,$^)) --o-composition-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/everything.%.PTB.qzv: Output/Step52/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "Premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/everything.%.Detail_PTB.qzv: Output/Step52/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "Detail Premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/everything.%.Simple_PTB.qzv: Output/Step52/everything.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "Simple Premature" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step52/%/ancom.tsv Output/Step52/%/data.tsv: Output/Step52/%.qzv | build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D))
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D)) 1> $(@D).stdout 2> $(@D).stderr

Output/Step52/%.pdf: Python/step08.py Output/Step52/%/ancom.tsv Output/Step52/%/data.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step52-DADA2-gg: Output/Step52/everything.DADA2.gg.qza Output/Step52/everything.DADA2.gg.pseudocount.qza Output/Step52/everything.DADA2.gg.pseudocount.PTB.qzv Output/Step52/everything.DADA2.gg.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.gg.pseudocount.PTB.pdf Output/Step52/everything.DADA2.gg.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.gg.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.gg.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.gg.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.gg.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.gg.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.gg.Mouth.qza Output/Step52/everything.DADA2.gg.Mouth.pseudocount.qza Output/Step52/everything.DADA2.gg.Mouth.pseudocount.PTB.qzv Output/Step52/everything.DADA2.gg.Mouth.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.gg.Mouth.pseudocount.PTB.pdf Output/Step52/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.gg.Mouth.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.gg.Mouth.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.gg.Mouth.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.gg.Mouth.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.gg.Cervix.qza Output/Step52/everything.DADA2.gg.Cervix.pseudocount.qza Output/Step52/everything.DADA2.gg.Cervix.pseudocount.PTB.qzv Output/Step52/everything.DADA2.gg.Cervix.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.gg.Cervix.pseudocount.PTB.pdf Output/Step52/everything.DADA2.gg.Cervix.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.gg.Cervix.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.gg.Cervix.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.gg.Cervix.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.gg.Cervix.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.gg.Cervix.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.gg.Vagina.qza Output/Step52/everything.DADA2.gg.Vagina.pseudocount.qza Output/Step52/everything.DADA2.gg.Vagina.pseudocount.PTB.qzv Output/Step52/everything.DADA2.gg.Vagina.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.gg.Vagina.pseudocount.PTB.pdf Output/Step52/everything.DADA2.gg.Vagina.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.gg.Vagina.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.gg.Vagina.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.gg.Vagina.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.gg.Vagina.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.gg.Vagina.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.gg.Neonate-1day.qza Output/Step52/everything.DADA2.gg.Neonate-1day.pseudocount.qza Output/Step52/everything.DADA2.gg.Neonate-1day.pseudocount.PTB.qzv Output/Step52/everything.DADA2.gg.Neonate-1day.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.gg.Neonate-1day.pseudocount.PTB.pdf Output/Step52/everything.DADA2.gg.Neonate-1day.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.gg.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.gg.Neonate-1day.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.gg.Neonate-1day.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.gg.Neonate-1day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.gg.Neonate-1day.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.gg.Neonate-3day.qza Output/Step52/everything.DADA2.gg.Neonate-3day.pseudocount.qza Output/Step52/everything.DADA2.gg.Neonate-3day.pseudocount.PTB.qzv Output/Step52/everything.DADA2.gg.Neonate-3day.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.gg.Neonate-3day.pseudocount.PTB.pdf Output/Step52/everything.DADA2.gg.Neonate-3day.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.gg.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.gg.Neonate-3day.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.gg.Neonate-3day.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.gg.Neonate-3day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.gg.Neonate-3day.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.gg.Neonate-3day.qza Output/Step52/everything.DADA2.gg.Neonate-5day.pseudocount.qza Output/Step52/everything.DADA2.gg.Neonate-5day.pseudocount.PTB.qzv Output/Step52/everything.DADA2.gg.Neonate-5day.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.gg.Neonate-5day.pseudocount.PTB.pdf Output/Step52/everything.DADA2.gg.Neonate-5day.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.gg.Neonate-5day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.gg.Neonate-5day.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.gg.Neonate-5day.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.gg.Neonate-5day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.gg.Neonate-5day.pseudocount.Simple_PTB.pdf
step52-DADA2-silva: Output/Step52/everything.DADA2.silva.qza Output/Step52/everything.DADA2.silva.pseudocount.qza Output/Step52/everything.DADA2.silva.pseudocount.PTB.qzv Output/Step52/everything.DADA2.silva.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.silva.pseudocount.PTB.pdf Output/Step52/everything.DADA2.silva.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.silva.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.silva.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.silva.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.silva.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.silva.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.silva.Mouth.qza Output/Step52/everything.DADA2.silva.Mouth.pseudocount.qza Output/Step52/everything.DADA2.silva.Mouth.pseudocount.PTB.qzv Output/Step52/everything.DADA2.silva.Mouth.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.silva.Mouth.pseudocount.PTB.pdf Output/Step52/everything.DADA2.silva.Mouth.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.silva.Mouth.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.silva.Mouth.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.silva.Mouth.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.silva.Mouth.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.silva.Mouth.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.silva.Cervix.qza Output/Step52/everything.DADA2.silva.Cervix.pseudocount.qza Output/Step52/everything.DADA2.silva.Cervix.pseudocount.PTB.qzv Output/Step52/everything.DADA2.silva.Cervix.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.silva.Cervix.pseudocount.PTB.pdf Output/Step52/everything.DADA2.silva.Cervix.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.silva.Cervix.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.silva.Cervix.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.silva.Cervix.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.silva.Cervix.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.silva.Cervix.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.silva.Vagina.qza Output/Step52/everything.DADA2.silva.Vagina.pseudocount.qza Output/Step52/everything.DADA2.silva.Vagina.pseudocount.PTB.qzv Output/Step52/everything.DADA2.silva.Vagina.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.silva.Vagina.pseudocount.PTB.pdf Output/Step52/everything.DADA2.silva.Vagina.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.silva.Vagina.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.silva.Vagina.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.silva.Vagina.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.silva.Vagina.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.silva.Vagina.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.silva.Neonate-1day.qza Output/Step52/everything.DADA2.silva.Neonate-1day.pseudocount.qza Output/Step52/everything.DADA2.silva.Neonate-1day.pseudocount.PTB.qzv Output/Step52/everything.DADA2.silva.Neonate-1day.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.silva.Neonate-1day.pseudocount.PTB.pdf Output/Step52/everything.DADA2.silva.Neonate-1day.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.silva.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.silva.Neonate-1day.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.silva.Neonate-1day.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.silva.Neonate-1day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.silva.Neonate-1day.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.silva.Neonate-3day.qza Output/Step52/everything.DADA2.silva.Neonate-3day.pseudocount.qza Output/Step52/everything.DADA2.silva.Neonate-3day.pseudocount.PTB.qzv Output/Step52/everything.DADA2.silva.Neonate-3day.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.silva.Neonate-3day.pseudocount.PTB.pdf Output/Step52/everything.DADA2.silva.Neonate-3day.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.silva.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.silva.Neonate-3day.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.silva.Neonate-3day.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.silva.Neonate-3day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.silva.Neonate-3day.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.silva.Neonate-3day.qza Output/Step52/everything.DADA2.silva.Neonate-5day.pseudocount.qza Output/Step52/everything.DADA2.silva.Neonate-5day.pseudocount.PTB.qzv Output/Step52/everything.DADA2.silva.Neonate-5day.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.silva.Neonate-5day.pseudocount.PTB.pdf Output/Step52/everything.DADA2.silva.Neonate-5day.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.silva.Neonate-5day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.silva.Neonate-5day.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.silva.Neonate-5day.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.silva.Neonate-5day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.silva.Neonate-5day.pseudocount.Simple_PTB.pdf
step52-DADA2-homd: Output/Step52/everything.DADA2.homd.qza Output/Step52/everything.DADA2.homd.pseudocount.qza Output/Step52/everything.DADA2.homd.pseudocount.PTB.qzv Output/Step52/everything.DADA2.homd.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.homd.pseudocount.PTB.pdf Output/Step52/everything.DADA2.homd.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.homd.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.homd.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.homd.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.homd.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.homd.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.homd.Mouth.qza Output/Step52/everything.DADA2.homd.Mouth.pseudocount.qza Output/Step52/everything.DADA2.homd.Mouth.pseudocount.PTB.qzv Output/Step52/everything.DADA2.homd.Mouth.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.homd.Mouth.pseudocount.PTB.pdf Output/Step52/everything.DADA2.homd.Mouth.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.homd.Mouth.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.homd.Mouth.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.homd.Mouth.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.homd.Mouth.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.homd.Mouth.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.homd.Cervix.qza Output/Step52/everything.DADA2.homd.Cervix.pseudocount.qza Output/Step52/everything.DADA2.homd.Cervix.pseudocount.PTB.qzv Output/Step52/everything.DADA2.homd.Cervix.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.homd.Cervix.pseudocount.PTB.pdf Output/Step52/everything.DADA2.homd.Cervix.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.homd.Cervix.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.homd.Cervix.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.homd.Cervix.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.homd.Cervix.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.homd.Cervix.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.homd.Vagina.qza Output/Step52/everything.DADA2.homd.Vagina.pseudocount.qza Output/Step52/everything.DADA2.homd.Vagina.pseudocount.PTB.qzv Output/Step52/everything.DADA2.homd.Vagina.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.homd.Vagina.pseudocount.PTB.pdf Output/Step52/everything.DADA2.homd.Vagina.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.homd.Vagina.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.homd.Vagina.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.homd.Vagina.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.homd.Vagina.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.homd.Vagina.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.homd.Neonate-1day.qza Output/Step52/everything.DADA2.homd.Neonate-1day.pseudocount.qza Output/Step52/everything.DADA2.homd.Neonate-1day.pseudocount.PTB.qzv Output/Step52/everything.DADA2.homd.Neonate-1day.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.homd.Neonate-1day.pseudocount.PTB.pdf Output/Step52/everything.DADA2.homd.Neonate-1day.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.homd.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.homd.Neonate-1day.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.homd.Neonate-1day.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.homd.Neonate-1day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.homd.Neonate-1day.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.homd.Neonate-3day.qza Output/Step52/everything.DADA2.homd.Neonate-3day.pseudocount.qza Output/Step52/everything.DADA2.homd.Neonate-3day.pseudocount.PTB.qzv Output/Step52/everything.DADA2.homd.Neonate-3day.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.homd.Neonate-3day.pseudocount.PTB.pdf Output/Step52/everything.DADA2.homd.Neonate-3day.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.homd.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.homd.Neonate-3day.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.homd.Neonate-3day.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.homd.Neonate-3day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.homd.Neonate-3day.pseudocount.Simple_PTB.pdf Output/Step52/everything.DADA2.homd.Neonate-3day.qza Output/Step52/everything.DADA2.homd.Neonate-5day.pseudocount.qza Output/Step52/everything.DADA2.homd.Neonate-5day.pseudocount.PTB.qzv Output/Step52/everything.DADA2.homd.Neonate-5day.pseudocount.PTB/data.tsv Output/Step52/everything.DADA2.homd.Neonate-5day.pseudocount.PTB.pdf Output/Step52/everything.DADA2.homd.Neonate-5day.pseudocount.Detail_PTB.qzv Output/Step52/everything.DADA2.homd.Neonate-5day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.DADA2.homd.Neonate-5day.pseudocount.Detail_PTB.pdf Output/Step52/everything.DADA2.homd.Neonate-5day.pseudocount.Simple_PTB.qzv Output/Step52/everything.DADA2.homd.Neonate-5day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.DADA2.homd.Neonate-5day.pseudocount.Simple_PTB.pdf
step52-Deblur-gg: Output/Step52/everything.Deblur.gg.qza Output/Step52/everything.Deblur.gg.pseudocount.qza Output/Step52/everything.Deblur.gg.pseudocount.PTB.qzv Output/Step52/everything.Deblur.gg.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.gg.pseudocount.PTB.pdf Output/Step52/everything.Deblur.gg.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.gg.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.gg.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.gg.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.gg.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.gg.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.gg.Mouth.qza Output/Step52/everything.Deblur.gg.Mouth.pseudocount.qza Output/Step52/everything.Deblur.gg.Mouth.pseudocount.PTB.qzv Output/Step52/everything.Deblur.gg.Mouth.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.gg.Mouth.pseudocount.PTB.pdf Output/Step52/everything.Deblur.gg.Mouth.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.gg.Mouth.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.gg.Mouth.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.gg.Mouth.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.gg.Mouth.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.gg.Mouth.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.gg.Cervix.qza Output/Step52/everything.Deblur.gg.Cervix.pseudocount.qza Output/Step52/everything.Deblur.gg.Cervix.pseudocount.PTB.qzv Output/Step52/everything.Deblur.gg.Cervix.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.gg.Cervix.pseudocount.PTB.pdf Output/Step52/everything.Deblur.gg.Cervix.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.gg.Cervix.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.gg.Cervix.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.gg.Cervix.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.gg.Cervix.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.gg.Cervix.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.gg.Vagina.qza Output/Step52/everything.Deblur.gg.Vagina.pseudocount.qza Output/Step52/everything.Deblur.gg.Vagina.pseudocount.PTB.qzv Output/Step52/everything.Deblur.gg.Vagina.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.gg.Vagina.pseudocount.PTB.pdf Output/Step52/everything.Deblur.gg.Vagina.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.gg.Vagina.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.gg.Vagina.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.gg.Vagina.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.gg.Vagina.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.gg.Vagina.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.gg.Neonate-1day.qza Output/Step52/everything.Deblur.gg.Neonate-1day.pseudocount.qza Output/Step52/everything.Deblur.gg.Neonate-1day.pseudocount.PTB.qzv Output/Step52/everything.Deblur.gg.Neonate-1day.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.gg.Neonate-1day.pseudocount.PTB.pdf Output/Step52/everything.Deblur.gg.Neonate-1day.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.gg.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.gg.Neonate-1day.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.gg.Neonate-1day.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.gg.Neonate-1day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.gg.Neonate-1day.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.gg.Neonate-3day.qza Output/Step52/everything.Deblur.gg.Neonate-3day.pseudocount.qza Output/Step52/everything.Deblur.gg.Neonate-3day.pseudocount.PTB.qzv Output/Step52/everything.Deblur.gg.Neonate-3day.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.gg.Neonate-3day.pseudocount.PTB.pdf Output/Step52/everything.Deblur.gg.Neonate-3day.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.gg.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.gg.Neonate-3day.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.gg.Neonate-3day.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.gg.Neonate-3day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.gg.Neonate-3day.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.gg.Neonate-5day.qza Output/Step52/everything.Deblur.gg.Neonate-5day.pseudocount.qza Output/Step52/everything.Deblur.gg.Neonate-5day.pseudocount.PTB.qzv Output/Step52/everything.Deblur.gg.Neonate-5day.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.gg.Neonate-5day.pseudocount.PTB.pdf Output/Step52/everything.Deblur.gg.Neonate-5day.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.gg.Neonate-5day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.gg.Neonate-5day.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.gg.Neonate-5day.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.gg.Neonate-5day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.gg.Neonate-5day.pseudocount.Simple_PTB.pdf
step52-Deblur-silva: Output/Step52/everything.Deblur.silva.qza Output/Step52/everything.Deblur.silva.pseudocount.qza Output/Step52/everything.Deblur.silva.pseudocount.PTB.qzv Output/Step52/everything.Deblur.silva.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.silva.pseudocount.PTB.pdf Output/Step52/everything.Deblur.silva.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.silva.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.silva.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.silva.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.silva.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.silva.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.silva.Mouth.qza Output/Step52/everything.Deblur.silva.Mouth.pseudocount.qza Output/Step52/everything.Deblur.silva.Mouth.pseudocount.PTB.qzv Output/Step52/everything.Deblur.silva.Mouth.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.silva.Mouth.pseudocount.PTB.pdf Output/Step52/everything.Deblur.silva.Mouth.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.silva.Mouth.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.silva.Mouth.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.silva.Mouth.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.silva.Mouth.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.silva.Mouth.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.silva.Cervix.qza Output/Step52/everything.Deblur.silva.Cervix.pseudocount.qza Output/Step52/everything.Deblur.silva.Cervix.pseudocount.PTB.qzv Output/Step52/everything.Deblur.silva.Cervix.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.silva.Cervix.pseudocount.PTB.pdf Output/Step52/everything.Deblur.silva.Cervix.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.silva.Cervix.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.silva.Cervix.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.silva.Cervix.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.silva.Cervix.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.silva.Cervix.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.silva.Vagina.qza Output/Step52/everything.Deblur.silva.Vagina.pseudocount.qza Output/Step52/everything.Deblur.silva.Vagina.pseudocount.PTB.qzv Output/Step52/everything.Deblur.silva.Vagina.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.silva.Vagina.pseudocount.PTB.pdf Output/Step52/everything.Deblur.silva.Vagina.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.silva.Vagina.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.silva.Vagina.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.silva.Vagina.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.silva.Vagina.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.silva.Vagina.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.silva.Neonate-1day.qza Output/Step52/everything.Deblur.silva.Neonate-1day.pseudocount.qza Output/Step52/everything.Deblur.silva.Neonate-1day.pseudocount.PTB.qzv Output/Step52/everything.Deblur.silva.Neonate-1day.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.silva.Neonate-1day.pseudocount.PTB.pdf Output/Step52/everything.Deblur.silva.Neonate-1day.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.silva.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.silva.Neonate-1day.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.silva.Neonate-1day.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.silva.Neonate-1day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.silva.Neonate-1day.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.silva.Neonate-3day.qza Output/Step52/everything.Deblur.silva.Neonate-3day.pseudocount.qza Output/Step52/everything.Deblur.silva.Neonate-3day.pseudocount.PTB.qzv Output/Step52/everything.Deblur.silva.Neonate-3day.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.silva.Neonate-3day.pseudocount.PTB.pdf Output/Step52/everything.Deblur.silva.Neonate-3day.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.silva.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.silva.Neonate-3day.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.silva.Neonate-3day.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.silva.Neonate-3day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.silva.Neonate-3day.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.silva.Neonate-3day.qza Output/Step52/everything.Deblur.silva.Neonate-5day.pseudocount.qza Output/Step52/everything.Deblur.silva.Neonate-5day.pseudocount.PTB.qzv Output/Step52/everything.Deblur.silva.Neonate-5day.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.silva.Neonate-5day.pseudocount.PTB.pdf Output/Step52/everything.Deblur.silva.Neonate-5day.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.silva.Neonate-5day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.silva.Neonate-5day.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.silva.Neonate-5day.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.silva.Neonate-5day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.silva.Neonate-5day.pseudocount.Simple_PTB.pdf
step52-Deblur-homd: Output/Step52/everything.Deblur.homd.qza Output/Step52/everything.Deblur.homd.pseudocount.qza Output/Step52/everything.Deblur.homd.pseudocount.PTB.qzv Output/Step52/everything.Deblur.homd.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.homd.pseudocount.PTB.pdf Output/Step52/everything.Deblur.homd.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.homd.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.homd.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.homd.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.homd.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.homd.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.homd.Mouth.qza Output/Step52/everything.Deblur.homd.Mouth.pseudocount.qza Output/Step52/everything.Deblur.homd.Mouth.pseudocount.PTB.qzv Output/Step52/everything.Deblur.homd.Mouth.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.homd.Mouth.pseudocount.PTB.pdf Output/Step52/everything.Deblur.homd.Mouth.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.homd.Mouth.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.homd.Mouth.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.homd.Mouth.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.homd.Mouth.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.homd.Mouth.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.homd.Cervix.qza Output/Step52/everything.Deblur.homd.Cervix.pseudocount.qza Output/Step52/everything.Deblur.homd.Cervix.pseudocount.PTB.qzv Output/Step52/everything.Deblur.homd.Cervix.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.homd.Cervix.pseudocount.PTB.pdf Output/Step52/everything.Deblur.homd.Cervix.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.homd.Cervix.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.homd.Cervix.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.homd.Cervix.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.homd.Cervix.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.homd.Cervix.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.homd.Vagina.qza Output/Step52/everything.Deblur.homd.Vagina.pseudocount.qza Output/Step52/everything.Deblur.homd.Vagina.pseudocount.PTB.qzv Output/Step52/everything.Deblur.homd.Vagina.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.homd.Vagina.pseudocount.PTB.pdf Output/Step52/everything.Deblur.homd.Vagina.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.homd.Vagina.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.homd.Vagina.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.homd.Vagina.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.homd.Vagina.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.homd.Vagina.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.homd.Neonate-1day.qza Output/Step52/everything.Deblur.homd.Neonate-1day.pseudocount.qza Output/Step52/everything.Deblur.homd.Neonate-1day.pseudocount.PTB.qzv Output/Step52/everything.Deblur.homd.Neonate-1day.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.homd.Neonate-1day.pseudocount.PTB.pdf Output/Step52/everything.Deblur.homd.Neonate-1day.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.homd.Neonate-1day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.homd.Neonate-1day.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.homd.Neonate-1day.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.homd.Neonate-1day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.homd.Neonate-1day.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.homd.Neonate-3day.qza Output/Step52/everything.Deblur.homd.Neonate-3day.pseudocount.qza Output/Step52/everything.Deblur.homd.Neonate-3day.pseudocount.PTB.qzv Output/Step52/everything.Deblur.homd.Neonate-3day.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.homd.Neonate-3day.pseudocount.PTB.pdf Output/Step52/everything.Deblur.homd.Neonate-3day.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.homd.Neonate-3day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.homd.Neonate-3day.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.homd.Neonate-3day.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.homd.Neonate-3day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.homd.Neonate-3day.pseudocount.Simple_PTB.pdf Output/Step52/everything.Deblur.homd.Neonate-3day.qza Output/Step52/everything.Deblur.homd.Neonate-5day.pseudocount.qza Output/Step52/everything.Deblur.homd.Neonate-5day.pseudocount.PTB.qzv Output/Step52/everything.Deblur.homd.Neonate-5day.pseudocount.PTB/data.tsv Output/Step52/everything.Deblur.homd.Neonate-5day.pseudocount.PTB.pdf Output/Step52/everything.Deblur.homd.Neonate-5day.pseudocount.Detail_PTB.qzv Output/Step52/everything.Deblur.homd.Neonate-5day.pseudocount.Detail_PTB/data.tsv Output/Step52/everything.Deblur.homd.Neonate-5day.pseudocount.Detail_PTB.pdf Output/Step52/everything.Deblur.homd.Neonate-5day.pseudocount.Simple_PTB.qzv Output/Step52/everything.Deblur.homd.Neonate-5day.pseudocount.Simple_PTB/data.tsv Output/Step52/everything.Deblur.homd.Neonate-5day.pseudocount.Simple_PTB.pdf
step52: step52-DADA2-gg step52-DADA2-silva step52-DADA2-homd step52-Deblur-gg step52-Deblur-silva step52-Deblur-homd
.PHONY: step52-DADA2-gg step52-DADA2-silva step52-DADA2-homd step52-Deblur-gg step52-Deblur-silva step52-Deblur-homd step52

# Step 53 (Abundance/Proportion Distribution)
Output/Step53 Output/Step53_Proportion:
	mkdir -p $@

Output/Step53/%.pdf: Python/step53.py Output/TSV/%.tsv Output/Metadata/everything.tsv | build.log Output/Step53
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step53_Proportion/%.pdf: Python/step53-2.py Output/TSV_Proportion/%.tsv Output/Metadata/everything.tsv | build.log Output/Step53_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step53-uncorrected: Output/Step53/everything.DADA2.gg.uncorrected.pdf Output/Step53/everything.DADA2.silva.uncorrected.pdf Output/Step53/everything.DADA2.homd.uncorrected.pdf Output/Step53/everything.Deblur.gg.uncorrected.pdf Output/Step53/everything.Deblur.silva.uncorrected.pdf Output/Step53/everything.Deblur.homd.uncorrected.pdf
step53-corrected: Output/Step53/everything.DADA2.gg.pdf Output/Step53/everything.DADA2.silva.pdf Output/Step53/everything.DADA2.homd.pdf Output/Step53/everything.Deblur.gg.pdf Output/Step53/everything.Deblur.silva.pdf Output/Step53/everything.Deblur.homd.pdf
step53-uncorrected-Proportion: Output/Step53_Proportion/everything.DADA2.gg.uncorrected.pdf Output/Step53_Proportion/everything.DADA2.silva.uncorrected.pdf Output/Step53_Proportion/everything.DADA2.homd.uncorrected.pdf Output/Step53_Proportion/everything.Deblur.gg.uncorrected.pdf Output/Step53_Proportion/everything.Deblur.silva.uncorrected.pdf Output/Step53_Proportion/everything.Deblur.homd.uncorrected.pdf
step53-corrected-Proportion: Output/Step53_Proportion/everything.DADA2.gg.pdf Output/Step53_Proportion/everything.DADA2.silva.pdf Output/Step53_Proportion/everything.DADA2.homd.pdf Output/Step53_Proportion/everything.Deblur.gg.pdf Output/Step53_Proportion/everything.Deblur.silva.pdf Output/Step53_Proportion/everything.Deblur.homd.pdf
step53: step53-uncorrected step53-corrected step53-uncorrected-Proportion step53-corrected-Proportion
.PHONY: step53-uncorrected step53-corrected step53-uncorrected-Proportion step53-corrected-Proportion step53

# Step 54 (LEfSe analysis after batch correction)
Output/Step54 Output/Step54/Default Output/Step54/Cladogram:
	mkdir -p $@
	chmod 777 $@

Output/Step54/%.DADA2.gg.Cervix.tsv Output/Step54/%.DADA2.gg.Mouth.tsv Output/Step54/%.DADA2.gg.Neonate-3day.tsv Output/Step54/%.DADA2.gg.Neonate-5day.tsv Output/Step54/%.DADA2.gg.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.gg.uncorrected.tsv Output/Metadata/%.tsv | Output/Step54 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step54/%.DADA2.silva.Cervix.tsv Output/Step54/%.DADA2.silva.Mouth.tsv Output/Step54/%.DADA2.silva.Neonate-3day.tsv Output/Step54/%.DADA2.silva.Neonate-5day.tsv Output/Step54/%.DADA2.silva.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.silva.uncorrected.tsv Output/Metadata/%.tsv | Output/Step54 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step54/%.DADA2.homd.Cervix.tsv Output/Step54/%.DADA2.homd.Mouth.tsv Output/Step54/%.DADA2.homd.Neonate-3day.tsv Output/Step54/%.DADA2.homd.Neonate-5day.tsv Output/Step54/%.DADA2.homd.Vagina.tsv: Python/step33.py Output/TSV/%.DADA2.homd.uncorrected.tsv Output/Metadata/%.tsv | Output/Step54 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step54/%.Deblur.gg.Cervix.tsv Output/Step54/%.Deblur.gg.Mouth.tsv Output/Step54/%.Deblur.gg.Neonate-3day.tsv Output/Step54/%.Deblur.gg.Neonate-5day.tsv Output/Step54/%.Deblur.gg.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.gg.uncorrected.tsv Output/Metadata/%.tsv | Output/Step54 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step54/%.Deblur.silva.Cervix.tsv Output/Step54/%.Deblur.silva.Mouth.tsv Output/Step54/%.Deblur.silva.Neonate-3day.tsv Output/Step54/%.Deblur.silva.Neonate-5day.tsv Output/Step54/%.Deblur.silva.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.silva.uncorrected.tsv Output/Metadata/%.tsv | Output/Step54 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step54/%.Deblur.homd.Cervix.tsv Output/Step54/%.Deblur.homd.Mouth.tsv Output/Step54/%.Deblur.homd.Neonate-3day.tsv Output/Step54/%.Deblur.homd.Neonate-5day.tsv Output/Step54/%.Deblur.homd.Vagina.tsv: Python/step33.py Output/TSV/%.Deblur.homd.uncorrected.tsv Output/Metadata/%.tsv | Output/Step54 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $(basename $@))) 1> $@.stdout 2> $@.stderr

Output/Step54/%.in: Output/Step54/%.tsv | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse format_input.py $(addprefix /,$^ $@) -c 2 -u 1 1> $@.stdout 2> $@.stderr

Output/Step54/%.res: Output/Step54/%.in | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse run_lefse.py $(addprefix /,$^ $@) --verbose 1 -y 1 1> $@.stdout 2> $@.stderr

Output/Step54/Default/%.pdf: Output/Step54/%.res | build.log Output/Step54/Default
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse plot_res.py $(addprefix /,$^ $@) --format pdf 1> $@.stdout 2> $@.stderr || $(DOCKER) touch $(addprefix /,$@)

Output/Step54/Cladogram/%.res: Output/Step54/%.res | build.log Output/Step54/Cladogram
	grep --perl-regexp "True|False" $^ > $@ || cp -f $^ $@

Output/Step54/Cladogram/%.pdf: Output/Step54/Cladogram/%.res | build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) biobakery/lefse plot_cladogram.py $(addprefix /,$^ $@) --format pdf 1> $@.stdout 2> $@.stderr || $(DOCKER) touch $(addprefix /,$@)

step54-tsv: Output/Step54/everything.DADA2.gg.Cervix.tsv Output/Step54/everything.DADA2.silva.Cervix.tsv Output/Step54/everything.DADA2.homd.Cervix.tsv Output/Step54/everything.Deblur.gg.Cervix.tsv Output/Step54/everything.Deblur.silva.Cervix.tsv
step54-in: $(addprefix Output/Step54/everything.DADA2.gg.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/everything.DADA2.silva.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/everything.DADA2.homd.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/everything.Deblur.gg.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/everything.Deblur.silva.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/everything.Deblur.homd.,$(addsuffix .in,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step54-res: $(addprefix Output/Step54/everything.DADA2.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/everything.DADA2.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/everything.DADA2.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/everything.Deblur.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/everything.Deblur.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/everything.Deblur.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step54-default: $(addprefix Output/Step54/Default/everything.DADA2.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Default/everything.DADA2.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Default/everything.DADA2.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Default/everything.Deblur.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Default/everything.Deblur.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Default/everything.Deblur.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step54-cladogram: $(addprefix Output/Step54/Cladogram/everything.DADA2.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Cladogram/everything.DADA2.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Cladogram/everything.DADA2.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Cladogram/everything.Deblur.gg.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Cladogram/everything.Deblur.silva.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Cladogram/everything.Deblur.homd.,$(addsuffix .res,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Cladogram/everything.DADA2.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Cladogram/everything.DADA2.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Cladogram/everything.DADA2.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Cladogram/everything.Deblur.gg.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Cladogram/everything.Deblur.silva.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina)) $(addprefix Output/Step54/Cladogram/everything.Deblur.homd.,$(addsuffix .pdf,Cervix Mouth Neonate-3day Neonate-5day Vagina))
step54: step54-tsv step54-in step54-res step54-default step54-cladogram
.PHONY: step54-tsv step54-in step54-res step54-default step54-cladogram step54

# Step 55 (Clustermap)
Output/Step55 Output/Step55_Proportion:
	mkdir -p $@

Output/Step55/%.pdf: Python/step55.py Output/TSV/%.tsv | build.log Output/Step55
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step55_Proportion/%.pdf: Python/step55-2.py Output/TSV_Proportion/%.tsv | build.log Output/Step55_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step55-uncorrected: Output/Step55/everything.DADA2.gg.uncorrected.pdf Output/Step55/everything.DADA2.silva.uncorrected.pdf Output/Step55/everything.DADA2.homd.uncorrected.pdf Output/Step55/everything.Deblur.gg.uncorrected.pdf Output/Step55/everything.Deblur.silva.uncorrected.pdf Output/Step55/everything.Deblur.homd.uncorrected.pdf
step55-corrected: Output/Step55/everything.DADA2.gg.pdf Output/Step55/everything.DADA2.silva.pdf Output/Step55/everything.DADA2.homd.pdf Output/Step55/everything.Deblur.gg.pdf Output/Step55/everything.Deblur.silva.pdf Output/Step55/everything.Deblur.homd.pdf
step55-uncorrected-Proportion: Output/Step55_Proportion/everything.DADA2.gg.uncorrected.pdf Output/Step55_Proportion/everything.DADA2.silva.uncorrected.pdf Output/Step55_Proportion/everything.DADA2.homd.uncorrected.pdf Output/Step55_Proportion/everything.Deblur.gg.uncorrected.pdf Output/Step55_Proportion/everything.Deblur.silva.uncorrected.pdf Output/Step55_Proportion/everything.Deblur.homd.uncorrected.pdf
step55-corrected-Proportion: Output/Step55_Proportion/everything.DADA2.gg.pdf Output/Step55_Proportion/everything.DADA2.silva.pdf Output/Step55_Proportion/everything.DADA2.homd.pdf Output/Step55_Proportion/everything.Deblur.gg.pdf Output/Step55_Proportion/everything.Deblur.silva.pdf Output/Step55_Proportion/everything.Deblur.homd.pdf
step55: step55-uncorrected step55-corrected step55-uncorrected-Proportion step55-corrected-Proportion
.PHONY: step55-uncorrected step55-corrected step55-uncorrected-Proportion step55-corrected-Proportion step55

# Step 56 (K-Nearest Neighbors Classifier)
Output/Step56 Output/Step56_Proportion:
	mkdir -p $@

Output/Step56/KNN.%.tar: Python/step56.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step56 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step56/KNN-two.%.tar: Python/step56-1.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step56 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step56_Proportion/KNN.%.tar: Python/step56.py Output/TSV_Proportion/everything.%.tsv Output/Metadata/everything.tsv | Output/Step56_Proportion build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step56_Proportion/KNN-two.%.tar: Python/step56-1.py Output/TSV_Proportion/everything.%.tsv Output/Metadata/everything.tsv | Output/Step56_Proportion build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step56-uncorrected: Output/Step56/KNN.DADA2.gg.uncorrected.tar Output/Step56/KNN.DADA2.silva.uncorrected.tar Output/Step56/KNN.DADA2.homd.uncorrected.tar Output/Step56/KNN.Deblur.gg.uncorrected.tar Output/Step56/KNN.Deblur.silva.uncorrected.tar Output/Step56/KNN.Deblur.homd.uncorrected.tar Output/Step56/KNN-two.DADA2.gg.uncorrected.tar Output/Step56/KNN-two.DADA2.silva.uncorrected.tar Output/Step56/KNN-two.DADA2.homd.uncorrected.tar Output/Step56/KNN-two.Deblur.gg.uncorrected.tar Output/Step56/KNN-two.Deblur.silva.uncorrected.tar Output/Step56/KNN-two.Deblur.homd.uncorrected.tar
step56-corrected: Output/Step56/KNN.DADA2.gg.tar Output/Step56/KNN.DADA2.silva.tar Output/Step56/KNN.DADA2.homd.tar Output/Step56/KNN.Deblur.gg.tar Output/Step56/KNN.Deblur.silva.tar Output/Step56/KNN.Deblur.homd.tar Output/Step56/KNN-two.DADA2.gg.tar Output/Step56/KNN-two.DADA2.silva.tar Output/Step56/KNN-two.DADA2.homd.tar Output/Step56/KNN-two.Deblur.gg.tar Output/Step56/KNN-two.Deblur.silva.tar Output/Step56/KNN-two.Deblur.homd.tar
step56-uncorrected-Proportion: Output/Step56_Proportion/KNN.DADA2.gg.uncorrected.tar Output/Step56_Proportion/KNN.DADA2.silva.uncorrected.tar Output/Step56_Proportion/KNN.DADA2.homd.uncorrected.tar Output/Step56_Proportion/KNN.Deblur.gg.uncorrected.tar Output/Step56_Proportion/KNN.Deblur.silva.uncorrected.tar Output/Step56_Proportion/KNN.Deblur.homd.uncorrected.tar Output/Step56_Proportion/KNN-two.DADA2.gg.uncorrected.tar Output/Step56_Proportion/KNN-two.DADA2.silva.uncorrected.tar Output/Step56_Proportion/KNN-two.DADA2.homd.uncorrected.tar Output/Step56_Proportion/KNN-two.Deblur.gg.uncorrected.tar Output/Step56_Proportion/KNN-two.Deblur.silva.uncorrected.tar Output/Step56_Proportion/KNN-two.Deblur.homd.uncorrected.tar
step56-corrected-Proportion: Output/Step56_Proportion/KNN.DADA2.gg.tar Output/Step56_Proportion/KNN.DADA2.silva.tar Output/Step56_Proportion/KNN.DADA2.homd.tar Output/Step56_Proportion/KNN.Deblur.gg.tar Output/Step56_Proportion/KNN.Deblur.silva.tar Output/Step56_Proportion/KNN.Deblur.homd.tar Output/Step56_Proportion/KNN-two.DADA2.gg.tar Output/Step56_Proportion/KNN-two.DADA2.silva.tar Output/Step56_Proportion/KNN-two.DADA2.homd.tar Output/Step56_Proportion/KNN-two.Deblur.gg.tar Output/Step56_Proportion/KNN-two.Deblur.silva.tar Output/Step56_Proportion/KNN-two.Deblur.homd.tar
step56: step56-uncorrected step56-corrected step56-uncorrected-Proportion step56-corrected-Proportion
.PHONY: step56-uncorrected step56-corrected step56-uncorrected-Proportion step56-corrected-Proportion step56

# Step 57 (Support Vector Machine Classifier with Bagging)
Output/Step57 Output/Step57_Proportion:
	mkdir -p $@

Output/Step57/SVM.%.tar: Python/step57.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step57 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step57/SVM-two.%.tar: Python/step57-1.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step57 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step57_Proportion/SVM.%.tar: Python/step57.py Output/TSV_Proportion/everything.%.tsv Output/Metadata/everything.tsv | Output/Step57_Proportion build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step57_Proportion/SVM-two.%.tar: Python/step57-1.py Output/TSV_Proportion/everything.%.tsv Output/Metadata/everything.tsv | Output/Step57_Proportion build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step57-uncorrected: Output/Step57/SVM.DADA2.gg.uncorrected.tar Output/Step57/SVM.DADA2.silva.uncorrected.tar Output/Step57/SVM.DADA2.homd.uncorrected.tar Output/Step57/SVM.Deblur.gg.uncorrected.tar Output/Step57/SVM.Deblur.silva.uncorrected.tar Output/Step57/SVM.Deblur.homd.uncorrected.tar Output/Step57/SVM-two.DADA2.gg.uncorrected.tar Output/Step57/SVM-two.DADA2.silva.uncorrected.tar Output/Step57/SVM-two.DADA2.homd.uncorrected.tar Output/Step57/SVM-two.Deblur.gg.uncorrected.tar Output/Step57/SVM-two.Deblur.silva.uncorrected.tar Output/Step57/SVM-two.Deblur.homd.uncorrected.tar
step57-corrected: Output/Step57/SVM.DADA2.gg.tar Output/Step57/SVM.DADA2.silva.tar Output/Step57/SVM.DADA2.homd.tar Output/Step57/SVM.Deblur.gg.tar Output/Step57/SVM.Deblur.silva.tar Output/Step57/SVM.Deblur.homd.tar Output/Step57/SVM-two.DADA2.gg.tar Output/Step57/SVM-two.DADA2.silva.tar Output/Step57/SVM-two.DADA2.homd.tar Output/Step57/SVM-two.Deblur.gg.tar Output/Step57/SVM-two.Deblur.silva.tar Output/Step57/SVM-two.Deblur.homd.tar
step57-uncorrected-Proportion: Output/Step57_Proportion/SVM.DADA2.gg.uncorrected.tar Output/Step57_Proportion/SVM.DADA2.silva.uncorrected.tar Output/Step57_Proportion/SVM.DADA2.homd.uncorrected.tar Output/Step57_Proportion/SVM.Deblur.gg.uncorrected.tar Output/Step57_Proportion/SVM.Deblur.silva.uncorrected.tar Output/Step57_Proportion/SVM.Deblur.homd.uncorrected.tar Output/Step57_Proportion/SVM-two.DADA2.gg.uncorrected.tar Output/Step57_Proportion/SVM-two.DADA2.silva.uncorrected.tar Output/Step57_Proportion/SVM-two.DADA2.homd.uncorrected.tar Output/Step57_Proportion/SVM-two.Deblur.gg.uncorrected.tar Output/Step57_Proportion/SVM-two.Deblur.silva.uncorrected.tar Output/Step57_Proportion/SVM-two.Deblur.homd.uncorrected.tar
step57-corrected-Proportion: Output/Step57_Proportion/SVM.DADA2.gg.tar Output/Step57_Proportion/SVM.DADA2.silva.tar Output/Step57_Proportion/SVM.DADA2.homd.tar Output/Step57_Proportion/SVM.Deblur.gg.tar Output/Step57_Proportion/SVM.Deblur.silva.tar Output/Step57_Proportion/SVM.Deblur.homd.tar Output/Step57_Proportion/SVM-two.DADA2.gg.tar Output/Step57_Proportion/SVM-two.DADA2.silva.tar Output/Step57_Proportion/SVM-two.DADA2.homd.tar Output/Step57_Proportion/SVM-two.Deblur.gg.tar Output/Step57_Proportion/SVM-two.Deblur.silva.tar Output/Step57_Proportion/SVM-two.Deblur.homd.tar
step57: step57-uncorrected step57-corrected step57-uncorrected-Proportion step57-corrected-Proportion
.PHONY: step57-uncorrected step57-corrected step57-uncorrected-Proportion step57-corrected-Proportion step57

# Step 58 (RandomForest Classifier with Bagging)
Output/Step58 Output/Step58_Proportion:
	mkdir -p $@

Output/Step58/RF.%.tar: Python/step58.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step58 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step58/RF-two.%.tar: Python/step58-1.py Output/TSV/everything.%.tsv Output/Metadata/everything.tsv | Output/Step58 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step58_Proportion/RF.%.tar: Python/step58.py Output/TSV_Proportion/everything.%.tsv Output/Metadata/everything.tsv | Output/Step58_Proportion build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step58_Proportion/RF-two.%.tar: Python/step58-1.py Output/TSV_Proportion/everything.%.tsv Output/Metadata/everything.tsv | Output/Step58_Proportion build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step58-uncorrected: Output/Step58/RF.DADA2.gg.uncorrected.tar Output/Step58/RF.DADA2.silva.uncorrected.tar Output/Step58/RF.DADA2.homd.uncorrected.tar Output/Step58/RF.Deblur.gg.uncorrected.tar Output/Step58/RF.Deblur.silva.uncorrected.tar Output/Step58/RF.Deblur.homd.uncorrected.tar Output/Step58/RF-two.DADA2.gg.uncorrected.tar Output/Step58/RF-two.DADA2.silva.uncorrected.tar Output/Step58/RF-two.DADA2.homd.uncorrected.tar Output/Step58/RF-two.Deblur.gg.uncorrected.tar Output/Step58/RF-two.Deblur.silva.uncorrected.tar Output/Step58/RF-two.Deblur.homd.uncorrected.tar
step58-corrected: Output/Step58/RF.DADA2.gg.tar Output/Step58/RF.DADA2.silva.tar Output/Step58/RF.DADA2.homd.tar Output/Step58/RF.Deblur.gg.tar Output/Step58/RF.Deblur.silva.tar Output/Step58/RF.Deblur.homd.tar Output/Step58/RF-two.DADA2.gg.tar Output/Step58/RF-two.DADA2.silva.tar Output/Step58/RF-two.DADA2.homd.tar Output/Step58/RF-two.Deblur.gg.tar Output/Step58/RF-two.Deblur.silva.tar Output/Step58/RF-two.Deblur.homd.tar
step58-uncorrected-Proportion: Output/Step58_Proportion/RF.DADA2.gg.uncorrected.tar Output/Step58_Proportion/RF.DADA2.silva.uncorrected.tar Output/Step58_Proportion/RF.DADA2.homd.uncorrected.tar Output/Step58_Proportion/RF.Deblur.gg.uncorrected.tar Output/Step58_Proportion/RF.Deblur.silva.uncorrected.tar Output/Step58_Proportion/RF.Deblur.homd.uncorrected.tar Output/Step58_Proportion/RF-two.DADA2.gg.uncorrected.tar Output/Step58_Proportion/RF-two.DADA2.silva.uncorrected.tar Output/Step58_Proportion/RF-two.DADA2.homd.uncorrected.tar Output/Step58_Proportion/RF-two.Deblur.gg.uncorrected.tar Output/Step58_Proportion/RF-two.Deblur.silva.uncorrected.tar Output/Step58_Proportion/RF-two.Deblur.homd.uncorrected.tar
step58-corrected-Proportion: Output/Step58_Proportion/RF.DADA2.gg.tar Output/Step58_Proportion/RF.DADA2.silva.tar Output/Step58_Proportion/RF.DADA2.homd.tar Output/Step58_Proportion/RF.Deblur.gg.tar Output/Step58_Proportion/RF.Deblur.silva.tar Output/Step58_Proportion/RF.Deblur.homd.tar Output/Step58_Proportion/RF-two.DADA2.gg.tar Output/Step58_Proportion/RF-two.DADA2.silva.tar Output/Step58_Proportion/RF-two.DADA2.homd.tar Output/Step58_Proportion/RF-two.Deblur.gg.tar Output/Step58_Proportion/RF-two.Deblur.silva.tar Output/Step58_Proportion/RF-two.Deblur.homd.tar
step58: step58-uncorrected step58-corrected step58-uncorrected-Proportion step58-corrected-Proportion
.PHONY: step58-uncorrected step58-corrected step58-uncorrected-Proportion step58-corrected-Proportion step58

# Step 59 (Violin Plots for Taxonomy Distribution)
Output/Step59 Output/Step59_Proportion:
	mkdir -p $@

Output/Step59/%.tar: Python/step59.py Output/TSV/%.tsv Output/Metadata/everything.tsv | Output/Step59 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step59_Proportion/%.tar: Python/step59.py Output/TSV_Proportion/%.tsv Output/Metadata/everything.tsv | Output/Step59_Proportion build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step59-uncorrected: Output/Step59/everything.DADA2.gg.uncorrected.tar Output/Step59/everything.DADA2.silva.uncorrected.tar Output/Step59/everything.DADA2.homd.uncorrected.tar Output/Step59/everything.Deblur.gg.uncorrected.tar Output/Step59/everything.Deblur.silva.uncorrected.tar Output/Step59/everything.Deblur.homd.uncorrected.tar
step59-corrected: Output/Step59/everything.DADA2.gg.tar Output/Step59/everything.DADA2.silva.tar Output/Step59/everything.DADA2.homd.tar Output/Step59/everything.Deblur.gg.tar Output/Step59/everything.Deblur.silva.tar Output/Step59/everything.Deblur.homd.tar
step59-uncorrected-Proportion: Output/Step59_Proportion/everything.DADA2.gg.uncorrected.tar Output/Step59_Proportion/everything.DADA2.silva.uncorrected.tar Output/Step59_Proportion/everything.DADA2.homd.uncorrected.tar Output/Step59_Proportion/everything.Deblur.gg.uncorrected.tar Output/Step59_Proportion/everything.Deblur.silva.uncorrected.tar Output/Step59_Proportion/everything.Deblur.homd.uncorrected.tar
step59-corrected-Proportion: Output/Step59_Proportion/everything.DADA2.gg.tar Output/Step59_Proportion/everything.DADA2.silva.tar Output/Step59_Proportion/everything.DADA2.homd.tar Output/Step59_Proportion/everything.Deblur.gg.tar Output/Step59_Proportion/everything.Deblur.silva.tar Output/Step59_Proportion/everything.Deblur.homd.tar
.PHONY: step59

# Step 60 (Abundance/Proportion Distribution)
Output/Step60 Output/Step60_Proportion:
	mkdir -p $@

Output/Step60/%.tar: Python/step60.py Output/TSV/%.tsv Output/Metadata/everything.tsv | Output/Step60 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step60_Proportion/%.tar: Python/step60.py Output/TSV_Proportion/%.tsv Output/Metadata/everything.tsv | Output/Step60_Proportion build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step60-uncorrected: Output/Step60/everything.DADA2.gg.uncorrected.tar Output/Step60/everything.DADA2.silva.uncorrected.tar Output/Step60/everything.DADA2.homd.uncorrected.tar Output/Step60/everything.Deblur.gg.uncorrected.tar Output/Step60/everything.Deblur.silva.uncorrected.tar Output/Step60/everything.Deblur.homd.uncorrected.tar
step60-corrected: Output/Step60/everything.DADA2.gg.tar Output/Step60/everything.DADA2.silva.tar Output/Step60/everything.DADA2.homd.tar Output/Step60/everything.Deblur.gg.tar Output/Step60/everything.Deblur.silva.tar Output/Step60/everything.Deblur.homd.tar
step60-uncorrected-Proportion: Output/Step60_Proportion/everything.DADA2.gg.uncorrected.tar Output/Step60_Proportion/everything.DADA2.silva.uncorrected.tar Output/Step60_Proportion/everything.DADA2.homd.uncorrected.tar Output/Step60_Proportion/everything.Deblur.gg.uncorrected.tar Output/Step60_Proportion/everything.Deblur.silva.uncorrected.tar Output/Step60_Proportion/everything.Deblur.homd.uncorrected.tar
step60-corrected-Proportion: Output/Step60_Proportion/everything.DADA2.gg.tar Output/Step60_Proportion/everything.DADA2.silva.tar Output/Step60_Proportion/everything.DADA2.homd.tar Output/Step60_Proportion/everything.Deblur.gg.tar Output/Step60_Proportion/everything.Deblur.silva.tar Output/Step60_Proportion/everything.Deblur.homd.tar
step60: step60-uncorrected step60-corrected step60-uncorrected-Proportion step60-corrected-Proportion
.PHONY: step60-uncorrected step60-corrected step60-uncorrected-Proportion step60-corrected-Proportion step60

# Step 61 (Violin plot with differentially abundant taxa from DESeq2)
Output/Step61:
	mkdir -p $@

Output/Step61/everything.%.tar: Python/step61.py Output/Step44/everything.%.tsv Output/Step44/everything.%.DAT.tsv Output/Step44/everything.%.coldata | build.log Output/Step61
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step61-uncorrected: $(addprefix Output/Step61/everything.DADA2.gg.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step61/everything.DADA2.silva.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step61/everything.DADA2.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step61/everything.Deblur.gg.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step61/everything.Deblur.silva.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step61/everything.Deblur.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina))
step61-corrected: $(addprefix Output/Step61/everything.DADA2.gg.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step61/everything.DADA2.silva.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step61/everything.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step61/everything.Deblur.gg.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step61/everything.Deblur.silva.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step61/everything.Deblur.homd.,$(addsuffix .tar,Cervix Mouth Vagina))
step61: step61-uncorrected step61-corrected
.PHONY: step61-uncorrected step61-corrected step61

# Step 62 (RandomForest Classifier with DAT)
Output/Step62 Output/Step62_Proportion:
	mkdir -p $@

Output/Step62/RF.%.tar: Python/step62.py Output/Step44/everything.%.tsv Output/Step44/everything.%.DAT.tsv Output/Metadata/everything.tsv | build.log Output/Step62
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step62/RF-two.%.tar: Python/step62-1.py Output/Step44/everything.%.tsv Output/Step44/everything.%.DAT.tsv Output/Metadata/everything.tsv | build.log Output/Step62
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step62_Proportion/RF.%.tar: Python/step62.py Output/Step44_Proportion/everything.%.tsv Output/Step44/everything.%.DAT.tsv Output/Metadata/everything.tsv | build.log Output/Step62_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step62_Proportion/RF-two.%.tar: Python/step62-1.py Output/Step44_Proportion/everything.%.tsv Output/Step44/everything.%.DAT.tsv Output/Metadata/everything.tsv | build.log Output/Step62_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step62-uncorrected: $(addprefix Output/Step62/RF.DADA2.gg.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF.DADA2.silva.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF.DADA2.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF.Deblur.gg.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF.Deblur.silva.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF.Deblur.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF-two.DADA2.gg.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF-two.DADA2.silva.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF-two.DADA2.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF-two.Deblur.gg.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF-two.Deblur.silva.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF-two.Deblur.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina))
step62-corrected: $(addprefix Output/Step62/RF.DADA2.gg.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF.DADA2.silva.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF.Deblur.gg.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF.Deblur.silva.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF.Deblur.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF-two.DADA2.gg.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF-two.DADA2.silva.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF-two.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF-two.Deblur.gg.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF-two.Deblur.silva.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62/RF-two.Deblur.homd.,$(addsuffix .tar,Cervix Mouth Vagina))
step62-uncorrected-proportion: $(addprefix Output/Step62_Proportion/RF.DADA2.gg.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF.DADA2.silva.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF.DADA2.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF.Deblur.gg.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF.Deblur.silva.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF.Deblur.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF-two.DADA2.gg.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF-two.DADA2.silva.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF-two.DADA2.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF-two.Deblur.gg.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF-two.Deblur.silva.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF-two.Deblur.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina))
step62-corrected-proportion: $(addprefix Output/Step62_Proportion/RF.DADA2.gg.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF.DADA2.silva.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF.Deblur.gg.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF.Deblur.silva.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF.Deblur.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF-two.DADA2.gg.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF-two.DADA2.silva.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF-two.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF-two.Deblur.gg.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF-two.Deblur.silva.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step62_Proportion/RF-two.Deblur.homd.,$(addsuffix .tar,Cervix Mouth Vagina))
step62: step62-uncorrected step62-corrected step62-uncorrected-proportion step62-corrected-proportion
.PHONY: step62-uncorrected step62-corrected step62-uncorrected-proportion step62-corrected-proportion step62

# Step 63 (Differentially abundant taxa with metagenomeSeq)
Output/Step63:
	mkdir -p $@

Output/Step63/everything.%.tsv: Python/step63-2.py Output/TSV/everything.%.tsv | build.log Output/Step63
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step63/metadata.%.tsv: Python/step63-1.py Output/Metadata/everything.tsv | build.log Output/Step63
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --site "$*" 1> $@.stdout 2> $@.stderr

Output/Step63/metadata.%.EL.tsv: Python/step63-1.py Output/Metadata/everything.tsv | build.log Output/Step63
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --site "$*" --EL 1> $@.stdout 2> $@.stderr

Output/Step63/metadata.%.EF.tsv: Python/step63-1.py Output/Metadata/everything.tsv | build.log Output/Step63
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --site "$*" --EF 1> $@.stdout 2> $@.stderr

Output/Step63/metadata.%.LF.tsv: Python/step63-1.py Output/Metadata/everything.tsv | build.log Output/Step63
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --site "$*" --LF 1> $@.stdout 2> $@.stderr

Output/Step63/metadata.%.PTB.tsv: Python/step63-1.py Output/Metadata/everything.tsv | build.log Output/Step63
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --site "$*" --PTB 1> $@.stdout 2> $@.stderr

Output/Step63/%.Mouth.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Mouth.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Cervix.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Cervix.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Vagina.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Vagina.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Mouth.EL.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Mouth.EL.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Cervix.EL.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Cervix.EL.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Vagina.EL.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Vagina.EL.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Mouth.EF.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Mouth.EF.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Cervix.EF.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Cervix.EF.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Vagina.EF.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Vagina.EF.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Mouth.LF.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Mouth.LF.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Cervix.LF.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Cervix.LF.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Vagina.LF.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Vagina.LF.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Mouth.PTB.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Mouth.PTB.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Cervix.PTB.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Cervix.PTB.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.Vagina.PTB.tsv: R/run_metagenomeSeq.R Output/Step63/%.uncorrected.tsv Output/Step63/metadata.Vagina.PTB.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(word 1,$(filter %.tsv,$^))) --metadata $(addprefix /,$(word 2,$(filter %.tsv,$^))) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step63/%.pdf: Python/step63.py Output/Step63/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step63-biom: Output/Step63/everything.DADA2.gg.uncorrected.tsv Output/Step63/everything.DADA2.silva.uncorrected.tsv Output/Step63/everything.DADA2.homd.uncorrected.tsv Output/Step63/everything.Deblur.gg.uncorrected.tsv Output/Step63/everything.Deblur.silva.uncorrected.tsv Output/Step63/everything.Deblur.homd.uncorrected.tsv
step63-metadata: Output/Step63/metadata.Mouth.tsv Output/Step63/metadata.Cervix.tsv Output/Step63/metadata.Vagina.tsv Output/Step63/metadata.Mouth.EL.tsv Output/Step63/metadata.Cervix.EL.tsv Output/Step63/metadata.Vagina.EL.tsv Output/Step63/metadata.Mouth.EF.tsv Output/Step63/metadata.Cervix.EF.tsv Output/Step63/metadata.Vagina.EF.tsv Output/Step63/metadata.Mouth.LF.tsv Output/Step63/metadata.Cervix.LF.tsv Output/Step63/metadata.Vagina.LF.tsv Output/Step63/metadata.Mouth.PTB.tsv Output/Step63/metadata.Cervix.PTB.tsv Output/Step63/metadata.Vagina.PTB.tsv
step63-tsv: Output/Step63/everything.DADA2.homd.Mouth.tsv Output/Step63/everything.DADA2.homd.Cervix.tsv Output/Step63/everything.DADA2.homd.Vagina.tsv Output/Step63/everything.DADA2.homd.Mouth.EL.tsv Output/Step63/everything.DADA2.homd.Cervix.EL.tsv Output/Step63/everything.DADA2.homd.Vagina.EL.tsv Output/Step63/everything.DADA2.homd.Mouth.EF.tsv Output/Step63/everything.DADA2.homd.Cervix.EF.tsv Output/Step63/everything.DADA2.homd.Vagina.EF.tsv Output/Step63/everything.DADA2.homd.Mouth.LF.tsv Output/Step63/everything.DADA2.homd.Cervix.LF.tsv Output/Step63/everything.DADA2.homd.Vagina.LF.tsv Output/Step63/everything.DADA2.homd.Mouth.PTB.tsv Output/Step63/everything.DADA2.homd.Cervix.PTB.tsv Output/Step63/everything.DADA2.homd.Vagina.PTB.tsv
step63-pdf: Output/Step63/everything.DADA2.homd.Mouth.pdf Output/Step63/everything.DADA2.homd.Cervix.pdf Output/Step63/everything.DADA2.homd.Vagina.pdf Output/Step63/everything.DADA2.homd.Mouth.EL.pdf Output/Step63/everything.DADA2.homd.Cervix.EL.pdf Output/Step63/everything.DADA2.homd.Vagina.EL.pdf Output/Step63/everything.DADA2.homd.Mouth.EF.pdf Output/Step63/everything.DADA2.homd.Cervix.EF.pdf Output/Step63/everything.DADA2.homd.Vagina.EF.pdf Output/Step63/everything.DADA2.homd.Mouth.LF.pdf Output/Step63/everything.DADA2.homd.Cervix.LF.pdf Output/Step63/everything.DADA2.homd.Vagina.LF.pdf Output/Step63/everything.DADA2.homd.Mouth.PTB.pdf Output/Step63/everything.DADA2.homd.Cervix.PTB.pdf Output/Step63/everything.DADA2.homd.Vagina.PTB.pdf
step63: step63-biom step63-metadata step63-tsv step63-pdf
.PHONY: step63-biom step63-metadata step63-tsv step63-pdf step63

# Step 64 (Differential abundance analysis with gneiss)
Output/Step64:
	mkdir -p $@

Output/Step64/everything.DADA2.PTB.qza: Output/DADA2/everything.table.qza | Output/Step64 build.log
	$(DOCKER) qiime gneiss correlation-clustering --i-table $(addprefix /,$(filter %.qza,$^)) --p-pseudocount 1.0 --o-clustering $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step64/everything.Deblur.PTB.qza: Output/Deblur/everything.table.qza | Output/Step64 build.log
	$(DOCKER) qiime gneiss correlation-clustering --i-table $(addprefix /,$(filter %.qza,$^)) --p-pseudocount 1.0 --o-clustering $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step64/everything.DADA2.GW.qza: Output/DADA2/everything.table.qza Output/Metadata/everything.tsv | Output/Step64 build.log
	$(DOCKER) qiime gneiss gradient-clustering --i-table $(addprefix /,$(filter %.qza,$^)) --m-gradient-file $(addprefix /,$(filter %.tsv,$^)) --m-gradient-column "Gestational Week" --p-ignore-missing-samples --o-clustering $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step64/everything.Deblur.GW.qza: Output/Deblur/everything.table.qza Output/Metadata/everything.tsv | Output/Step64 build.log
	$(DOCKER) qiime gneiss gradient-clustering --i-table $(addprefix /,$(filter %.qza,$^)) --m-gradient-file $(addprefix /,$(filter %.tsv,$^)) --m-gradient-column "Gestational Week" --p-ignore-missing-samples --o-clustering $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step64/homd.alignment.qza Output/Step64/homd.masked_alignment.qza Output/Step64/homd.tree.qza Output/Step64/homd.rooted_tree.qza: Output/Taxonomy/homd_fasta.qza | Output/Step64 build.log
	$(DOCKER) qiime phylogeny align-to-tree-mafft-fasttree --i-sequences $(addprefix /,$^) --p-n-threads $(CPUS) --o-alignment $(addprefix /,$(@D)/homd.alignment.qza) --o-masked-alignment $(addprefix /,$(@D)/homd.masked_alignment.qza) --o-tree $(addprefix /,$(@D)/homd.tree.qza) --o-rooted-tree $(addprefix /,$(@D)/homd.rooted_tree.qza) 1> $@.stdout 2> $@.stderr

Output/Step64/everything.DADA2.%.hierarchical.qza: Output/DADA2/everything.table.qza Output/Step64/everything.DADA2.%.qza | build.log
	$(DOCKER) qiime gneiss ilr-hierarchical --i-table $(addprefix /,$(word 1,$(filter %.qza,$^))) --i-tree $(addprefix /,$(word 2,$(filter %.qza,$^))) --p-pseudocount 1.0 --o-balances $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step64/everything.Deblur.%.hierarchical.qza: Output/Deblur/everything.table.qza Output/Step64/everything.Deblur.%.qza | build.log
	$(DOCKER) qiime gneiss ilr-hierarchical --i-table $(addprefix /,$(word 1,$(filter %.qza,$^))) --i-tree $(addprefix /,$(word 2,$(filter %.qza,$^))) --p-pseudocount 1.0 --o-balances $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step64/everything.DADA2.%.hierarchical.tsv: Bash/export_TSV.bash Output/Step64/everything.DADA2.%.hierarchical.qza Output/Taxonomy/everything.DADA2.homd.qza | build.log
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step64/everything.Deblur.%.hierarchical.tsv: Bash/export_TSV.bash Output/Step64/everything.Deblur.%.hierarchical.qza Output/Taxonomy/everything.Deblur.homd.qza | build.log
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step64/everything.DADA2.%.dendrogram.qzv: Output/DADA2/everything.table.qza Output/Step64/everything.DADA2.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime gneiss dendrogram-heatmap --i-table $(addprefix /,$(word 1,$(filter %.qza,$^))) --i-tree $(addprefix /,$(word 2,$(filter %.qza,$^))) --m-metadata-file $(addprefix /,$(filter %.tsv,$^)) --m-metadata-column "Detail Premature" --p-pseudocount 1.0 --p-method "log" --p-color-map "coolwarm" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Step64/everything.Deblur.%.dendrogram.qzv: Output/Deblur/everything.table.qza Output/Step64/everything.Deblur.%.qza Output/Metadata/everything.tsv | build.log
	$(DOCKER) qiime gneiss dendrogram-heatmap --i-table $(addprefix /,$(word 1,$(filter %.qza,$^))) --i-tree $(addprefix /,$(word 2,$(filter %.qza,$^))) --m-metadata-file $(addprefix /,$(filter %.tsv,$^)) --m-metadata-column "Detail Premature" --p-pseudocount 1.0 --p-method "log" --p-color-map "coolwarm" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step64-clustering: Output/Step64/everything.DADA2.PTB.qza Output/Step64/everything.Deblur.PTB.qza Output/Step64/everything.DADA2.GW.qza Output/Step64/everything.Deblur.GW.qza
step64-dendrogram: Output/Step64/everything.DADA2.PTB.dendrogram.qzv Output/Step64/everything.Deblur.PTB.dendrogram.qzv Output/Step64/everything.DADA2.GW.dendrogram.qzv Output/Step64/everything.Deblur.GW.dendrogram.qzv
step64: step64-clustering step64-dendrogram
.PHONY: step64-clustering step64-dendrogram step64

# Step 65 (Venn diagram for differentially abundant taxa of DESeq2)
Output/Step65 Output/Step65_Updown:
	mkdir -p $@

Output/Step65/everything.%.Mouth.Up.pdf: Python/step65.py Output/Step44/everything.%.Mouth.DAT.tsv Output/Step44/everything.%.EL.Mouth.DAT.tsv Output/Step44/everything.%.EF.Mouth.DAT.tsv Output/Step44/everything.%.LF.Mouth.DAT.tsv | build.log Output/Step65
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "Early PTB vs. Late PTB" "Earlry PTB vs. Normal" "Late PTB vs. Normal" --up 1> $@.stdout 2> $@.stderr

Output/Step65/everything.%.Mouth.Down.pdf: Python/step65.py Output/Step44/everything.%.Mouth.DAT.tsv Output/Step44/everything.%.EL.Mouth.DAT.tsv Output/Step44/everything.%.EF.Mouth.DAT.tsv Output/Step44/everything.%.LF.Mouth.DAT.tsv | build.log Output/Step65
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "Early PTB vs. Late PTB" "Earlry PTB vs. Normal" "Late PTB vs. Normal" --down 1> $@.stdout 2> $@.stderr

Output/Step65/everything.%.Cervix.Up.pdf: Python/step65.py Output/Step44/everything.%.Cervix.DAT.tsv Output/Step44/everything.%.EL.Cervix.DAT.tsv Output/Step44/everything.%.EF.Cervix.DAT.tsv Output/Step44/everything.%.LF.Cervix.DAT.tsv | build.log Output/Step65
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "Early PTB vs. Late PTB" "Earlry PTB vs. Normal" "Late PTB vs. Normal" --up 1> $@.stdout 2> $@.stderr

Output/Step65/everything.%.Cervix.Down.pdf: Python/step65.py Output/Step44/everything.%.Cervix.DAT.tsv Output/Step44/everything.%.EL.Cervix.DAT.tsv Output/Step44/everything.%.EF.Cervix.DAT.tsv Output/Step44/everything.%.LF.Cervix.DAT.tsv | build.log Output/Step65
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "Early PTB vs. Late PTB" "Earlry PTB vs. Normal" "Late PTB vs. Normal" --down 1> $@.stdout 2> $@.stderr

Output/Step65/everything.%.Vagina.Up.pdf: Python/step65.py Output/Step44/everything.%.Vagina.DAT.tsv Output/Step44/everything.%.EL.Vagina.DAT.tsv Output/Step44/everything.%.EF.Vagina.DAT.tsv Output/Step44/everything.%.LF.Vagina.DAT.tsv | build.log Output/Step65
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "Early PTB vs. Late PTB" "Earlry PTB vs. Normal" "Late PTB vs. Normal" --up 1> $@.stdout 2> $@.stderr

Output/Step65/everything.%.Vagina.Down.pdf: Python/step65.py Output/Step44/everything.%.Vagina.DAT.tsv Output/Step44/everything.%.EL.Vagina.DAT.tsv Output/Step44/everything.%.EF.Vagina.DAT.tsv Output/Step44/everything.%.LF.Vagina.DAT.tsv | build.log Output/Step65
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "Early PTB vs. Late PTB" "Earlry PTB vs. Normal" "Late PTB vs. Normal" --down 1> $@.stdout 2> $@.stderr

Output/Step65_Updown/everything.%.Mouth.tsv: Python/step65-1.py Output/Step44/everything.%.Mouth.DAT.tsv Output/Step44/everything.%.EL.Mouth.DAT.tsv Output/Step44/everything.%.EF.Mouth.DAT.tsv Output/Step44/everything.%.LF.Mouth.DAT.tsv | build.log Output/Step65_Updown
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "Early PTB vs. Late PTB" "Earlry PTB vs. Normal" "Late PTB vs. Normal" 1> $@.stdout 2> $@.stderr

Output/Step65_Updown/everything.%.Cervix.tsv: Python/step65-1.py Output/Step44/everything.%.Cervix.DAT.tsv Output/Step44/everything.%.EL.Cervix.DAT.tsv Output/Step44/everything.%.EF.Cervix.DAT.tsv Output/Step44/everything.%.LF.Cervix.DAT.tsv | build.log Output/Step65_Updown
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "Early PTB vs. Late PTB" "Earlry PTB vs. Normal" "Late PTB vs. Normal" 1> $@.stdout 2> $@.stderr

Output/Step65_Updown/everything.%.Vagina.tsv: Python/step65-1.py Output/Step44/everything.%.Vagina.DAT.tsv Output/Step44/everything.%.EL.Vagina.DAT.tsv Output/Step44/everything.%.EF.Vagina.DAT.tsv Output/Step44/everything.%.LF.Vagina.DAT.tsv | build.log Output/Step65_Updown
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "Early PTB vs. Late PTB" "Earlry PTB vs. Normal" "Late PTB vs. Normal" 1> $@.stdout 2> $@.stderr

step65-multiple: Output/Step65/everything.DADA2.homd.Mouth.Up.pdf Output/Step65/everything.DADA2.homd.Mouth.Down.pdf Output/Step65/everything.DADA2.homd.uncorrected.Mouth.Up.pdf Output/Step65/everything.DADA2.homd.uncorrected.Mouth.Down.pdf Output/Step65/everything.DADA2.homd.Cervix.Up.pdf Output/Step65/everything.DADA2.homd.Cervix.Down.pdf Output/Step65/everything.DADA2.homd.uncorrected.Cervix.Up.pdf Output/Step65/everything.DADA2.homd.uncorrected.Cervix.Down.pdf Output/Step65/everything.DADA2.homd.Vagina.Up.pdf Output/Step65/everything.DADA2.homd.Vagina.Down.pdf Output/Step65/everything.DADA2.homd.uncorrected.Vagina.Up.pdf Output/Step65/everything.DADA2.homd.uncorrected.Vagina.Down.pdf
step65-updown: Output/Step65_Updown/everything.DADA2.homd.Mouth.tsv Output/Step65_Updown/everything.DADA2.homd.uncorrected.Mouth.tsv Output/Step65_Updown/everything.DADA2.homd.Cervix.tsv Output/Step65_Updown/everything.DADA2.homd.uncorrected.Cervix.tsv Output/Step65_Updown/everything.DADA2.homd.Vagina.tsv Output/Step65_Updown/everything.DADA2.homd.uncorrected.Vagina.tsv
step65: step65-multiple step65-updown
.PHONY: step65-multiple step65-updown step65

# Step 66 (Select DAT of DESeq2 for classification)
Output/Step66 Output/Step66_Proportion:
	mkdir -p $@

Output/Step66/everything.%.Mouth.tsv: Python/step66.py Output/Step44/everything.%.Mouth.tsv Output/Step44/everything.%.Mouth.DAT.tsv Output/Step44/everything.%.EL.Mouth.DAT.tsv Output/Step44/everything.%.EF.Mouth.DAT.tsv Output/Step44/everything.%.LF.Mouth.DAT.tsv | build.log Output/Step66
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step66/everything.%.Vagina.tsv: Python/step66.py Output/Step44/everything.%.Vagina.tsv Output/Step44/everything.%.Vagina.DAT.tsv Output/Step44/everything.%.EL.Vagina.DAT.tsv Output/Step44/everything.%.EF.Vagina.DAT.tsv Output/Step44/everything.%.LF.Vagina.DAT.tsv | build.log Output/Step66
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step66/everything.%.Cervix.tsv: Python/step66.py Output/Step44/everything.%.Cervix.tsv Output/Step44/everything.%.Cervix.DAT.tsv Output/Step44/everything.%.EL.Cervix.DAT.tsv Output/Step44/everything.%.EF.Cervix.DAT.tsv Output/Step44/everything.%.LF.Cervix.DAT.tsv | build.log Output/Step66
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step66_Proportion/everything.%.Mouth.tsv: Python/step66.py Output/Step44_Proportion/everything.%.Mouth.tsv Output/Step44/everything.%.Mouth.DAT.tsv Output/Step44/everything.%.EL.Mouth.DAT.tsv Output/Step44/everything.%.EF.Mouth.DAT.tsv Output/Step44/everything.%.LF.Mouth.DAT.tsv | build.log Output/Step66_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step66_Proportion/everything.%.Vagina.tsv: Python/step66.py Output/Step44_Proportion/everything.%.Vagina.tsv Output/Step44/everything.%.Vagina.DAT.tsv Output/Step44/everything.%.EL.Vagina.DAT.tsv Output/Step44/everything.%.EF.Vagina.DAT.tsv Output/Step44/everything.%.LF.Vagina.DAT.tsv | build.log Output/Step66_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step66_Proportion/everything.%.Cervix.tsv: Python/step66.py Output/Step44_Proportion/everything.%.Cervix.tsv Output/Step44/everything.%.Cervix.DAT.tsv Output/Step44/everything.%.EL.Cervix.DAT.tsv Output/Step44/everything.%.EF.Cervix.DAT.tsv Output/Step44/everything.%.LF.Cervix.DAT.tsv | build.log Output/Step66_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step66: Output/Step66/everything.DADA2.homd.Mouth.tsv Output/Step66/everything.DADA2.homd.uncorrected.Mouth.tsv Output/Step66/everything.DADA2.homd.Vagina.tsv Output/Step66/everything.DADA2.homd.uncorrected.Vagina.tsv Output/Step66/everything.DADA2.homd.Cervix.tsv Output/Step66/everything.DADA2.homd.uncorrected.Cervix.tsv Output/Step66_Proportion/everything.DADA2.homd.Mouth.tsv Output/Step66_Proportion/everything.DADA2.homd.uncorrected.Mouth.tsv Output/Step66_Proportion/everything.DADA2.homd.Vagina.tsv Output/Step66_Proportion/everything.DADA2.homd.uncorrected.Vagina.tsv Output/Step66_Proportion/everything.DADA2.homd.Cervix.tsv Output/Step66_Proportion/everything.DADA2.homd.uncorrected.Cervix.tsv
.PHONY: step66

# Step 67 (RandomForest Classifier with DAT-union of DESeq2)
Output/Step67 Output/Step67_Proportion:
	mkdir -p $@

Output/Step67/RF.%.tar: Python/step67.py Output/Step66/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step67
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step67/RF-two.%.tar: Python/step67-1.py Output/Step66/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step67
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step67_Proportion/RF.%.tar: Python/step67.py Output/Step66_Proportion/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step67_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step67_Proportion/RF-two.%.tar: Python/step67-1.py Output/Step66_Proportion/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step67_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step67-uncorrected: $(addprefix Output/Step67/RF.DADA2.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step67/RF-two.DADA2.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina))
step67-corrected: $(addprefix Output/Step67/RF.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step67/RF-two.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina))
step67-uncorrected-proportion: $(addprefix Output/Step67_Proportion/RF.DADA2.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step67_Proportion/RF-two.DADA2.homd.uncorrected.,$(addsuffix .tar,Cervix Mouth Vagina))
step67-corrected-proportion: $(addprefix Output/Step67_Proportion/RF.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step67_Proportion/RF-two.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina))
step67: step67-uncorrected step67-corrected step67-uncorrected-proportion step67-corrected-proportion
.PHONY: step67-uncorrected step67-corrected step67-uncorrected-proportion step67-corrected-proportion step67

# Step 68 (Venn diagram for differentially abundant taxa of metagenomeSeq)
Output/Step68 Output/Step68_Updown:
	mkdir -p $@

Output/Step68/everything.%.Up.pdf: Python/step68.py Output/Step63/everything.%.tsv Output/Step63/everything.%.EL.tsv Output/Step63/everything.%.EF.tsv Output/Step63/everything.%.LF.tsv | build.log Output/Step68
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "Early PTB vs. Late PTB" "Early PTB vs. Normal" "Late PTB vs. Normal" --up 1> $@.stdout 2> $@.stderr

Output/Step68/everything.%.Down.pdf: Python/step68.py Output/Step63/everything.%.tsv Output/Step63/everything.%.EL.tsv Output/Step63/everything.%.EF.tsv Output/Step63/everything.%.LF.tsv | build.log Output/Step68
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "Early PTB vs. Late PTB" "Early PTB vs. Normal" "Late PTB vs. Normal" --down 1> $@.stdout 2> $@.stderr

Output/Step68_Updown/everything.%.pdf: Python/step68-1.py Output/Step63/everything.%.tsv Output/Step63/everything.%.PTB.tsv Output/Step63/everything.%.EL.tsv Output/Step63/everything.%.EF.tsv Output/Step63/everything.%.LF.tsv | build.log Output/Step68_Updown
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "PTB vs. Normal" "Early PTB vs. Late PTB" "Early PTB vs. Normal" "Late PTB vs. Normal" 1> $@.stdout 2> $@.stderr

Output/Step68_Updown/everything.%.tsv: Python/step68-2.py Output/Step63/everything.%.tsv Output/Step63/everything.%.PTB.tsv Output/Step63/everything.%.EL.tsv Output/Step63/everything.%.EF.tsv Output/Step63/everything.%.LF.tsv | build.log Output/Step68_Updown
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Early PTB vs. Late PTB+Normal" "PTB vs. Normal" "Early PTB vs. Late PTB" "Early PTB vs. Normal" "Late PTB vs. Normal" 1> $@.stdout 2> $@.stderr

step68-volcano: Output/Step68/everything.DADA2.homd.Mouth.Up.pdf Output/Step68/everything.DADA2.homd.Mouth.Down.pdf Output/Step68/everything.DADA2.homd.Cervix.Up.pdf Output/Step68/everything.DADA2.homd.Cervix.Down.pdf Output/Step68/everything.DADA2.homd.Vagina.Up.pdf Output/Step68/everything.DADA2.homd.Vagina.Down.pdf
step68-upset: Output/Step68_Updown/everything.DADA2.homd.Mouth.pdf Output/Step68_Updown/everything.DADA2.homd.Cervix.pdf Output/Step68_Updown/everything.DADA2.homd.Vagina.pdf
step68-tsv: Output/Step68_Updown/everything.DADA2.homd.Mouth.tsv Output/Step68_Updown/everything.DADA2.homd.Cervix.tsv Output/Step68_Updown/everything.DADA2.homd.Vagina.tsv
step68: step68-volcano step68-upset step68-tsv
.PHONY: step68-volcano step68-upset step68-tsv step68

# Step 69 (Select DAT of metagenomeSeq for classification)
Output/Step69 Output/Step69_Proportion:
	mkdir -p $@

Output/Step69/everything.%.tsv: Python/step69.py Output/Step44/everything.%.tsv Output/Step63/everything.%.tsv Output/Step63/everything.%.PTB.tsv Output/Step63/everything.%.EL.tsv Output/Step63/everything.%.EF.tsv Output/Step63/everything.%.LF.tsv | build.log Output/Step69
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step69_Proportion/everything.%.tsv: Python/step69.py Output/Step44_Proportion/everything.%.tsv Output/Step63/everything.%.tsv Output/Step63/everything.%.PTB.tsv Output/Step63/everything.%.EL.tsv Output/Step63/everything.%.EF.tsv Output/Step63/everything.%.LF.tsv | build.log Output/Step69_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step69: Output/Step69/everything.DADA2.homd.Mouth.tsv Output/Step69/everything.DADA2.homd.Vagina.tsv Output/Step69/everything.DADA2.homd.Cervix.tsv Output/Step69_Proportion/everything.DADA2.homd.Mouth.tsv Output/Step69_Proportion/everything.DADA2.homd.Vagina.tsv Output/Step69_Proportion/everything.DADA2.homd.Cervix.tsv
.PHONY: step69

# Step 70 (RandomForest Classifier with DAT-union of metagenomeSeq)
Output/Step70 Output/Step70_Proportion:
	mkdir -p $@

Output/Step70/RF.%.tar: Python/step67.py Output/Step69/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step70
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step70/RF-two.%.tar: Python/step67-1.py Output/Step69/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step70
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step70/RF-PTB.%.tar: Python/step67-2.py Output/Step69/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step70
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step70_Proportion/RF.%.tar: Python/step67.py Output/Step69_Proportion/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step70_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step70_Proportion/RF-two.%.tar: Python/step67-1.py Output/Step69_Proportion/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step70_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Step70_Proportion/RF-PTB.%.tar: Python/step67-2.py Output/Step69_Proportion/everything.%.tsv Output/Metadata/everything.tsv | build.log Output/Step70_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step70-corrected: $(addprefix Output/Step70/RF.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step70/RF-two.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step70/RF-PTB.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina))
step70-corrected-proportion: $(addprefix Output/Step70_Proportion/RF.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step70_Proportion/RF-two.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step70_Proportion/RF-PTB.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina))
step70: step70-corrected step70-corrected-proportion
.PHONY: step70-corrected step70-corrected-proportion step70

# Step 71 (Violin plot with differentially abundant taxa from metagenomeSeq)
Output/Step71 Output/Step71_Proportion:
	mkdir -p $@

Output/Step71/everything.%.tar: Python/step71.py Output/Step69/everything.%.tsv Output/Metadata/everything.tsv Output/Step63/everything.%.EL.tsv Output/Step63/everything.%.EF.tsv Output/Step63/everything.%.LF.tsv | build.log Output/Step71
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step71_Proportion/everything.%.tar: Python/step71.py Output/Step69_Proportion/everything.%.tsv Output/Metadata/everything.tsv Output/Step63/everything.%.EL.tsv Output/Step63/everything.%.EF.tsv Output/Step63/everything.%.LF.tsv | build.log Output/Step71_Proportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step71: $(addprefix Output/Step71/everything.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina)) $(addprefix Output/Step71_Proportion/everything.DADA2.homd.,$(addsuffix .tar,Cervix Mouth Vagina))
.PHONY: step71

# Step 72 (PICRUSt2)
Output/Step72:
	mkdir -p $@

Output/Step72/everything.DADA2.fna: Bash/export_fna.bash Output/DADA2/everything.seq.qza | build.log Output/Step72
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Step72/%.biom: Output/TSV/%.tsv | build.log Output/Step72
	$(DOCKER) biom convert --input-fp $(addprefix /,$^) --output-fp $(addprefix /,$@) --table-type="OTU table" --to-hdf5

Output/Step72/everything.DADA2.homd/path_abun_unstrat_per_seq.tsv.gz Output/Step72/everything.DADA2.homd/path_abun_unstrat.tsv.gz: Output/Step72/everything.DADA2.fna Output/Step72/everything.DADA2.homd.biom | build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D))
	$(DOCKER) picrust2_pipeline.py --study_fasta $(addprefix /,$(filter %.fna,$^)) --input $(addprefix /,$(filter %.biom,$^)) --output $(addprefix /,$(@D)) --processes $(CPUS) --stratified --per_sequence_contrib --remove_intermediate 1> $(@D).stdout 2> $(@D).stderr
	$(DOCKER) cp $(addprefix /,$(@D)/pathways_out/** $(@D))

Output/Step72/everything.DADA2.homd.%.tsv.gz: Output/Step72/everything.DADA2.homd/%.tsv.gz | build.log
	$(DOCKER) add_descriptions.py --input $(addprefix /,$^) --output $(addprefix /,$@) --map_type "METACYC" 1> $@.stdout 2> $@.stderr

step72: Output/Step72/everything.DADA2.fna Output/Step72/everything.DADA2.fna Output/Step72/everything.DADA2.homd.biom Output/Step72/everything.DADA2.homd/path_abun_unstrat_per_seq.tsv.gz Output/Step72/path_abun_unstrat_per_seq.tsv.gz Output/Step72/path_abun_unstrat.tsv.gz
.PHONY: step72
