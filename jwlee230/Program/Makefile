# Variables
DATE := $(shell date "+%y%m%d")
RANDOM := $(shell bash -c 'echo $$RANDOM')
IMAGE_NAME = premature:latest
CPUS = 10
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
latest: step08
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
	$(DOCKER) qiime dada2 denoise-paired --i-demultiplexed-seqs $(addprefix /,$<) --p-n-threads $(CPUS) --p-trunc-len-f 265 --p-trunc-len-r 265 --o-table $(addprefix /,$(basename $(basename $@))).table.qza --o-representative-sequences $(addprefix /,$(basename $(basename $@))).seq.qza --o-denoising-stats $(addprefix /,$(basename $(basename $@))).stat.qza 1> $@.stdout 2> $@.stderr

step04: Output/DADA2/whole.seq.qza
.PHONY += step04

# Step 05 (Deblur Denoising)
Output/Deblur:
	mkdir -p $@

Output/Deblur/%.table.qza Output/Deblur/%.seq.qza Output/Deblur/%.stat.qza: Output/Filter/%.seq.qza | Output/Deblur build.log
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

Output/Metadata/whole.tsv: Python/step07.py $(wildcard Data/RawData/*.fastq.gz) | Output/Metadata build.log
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

Output/ANCOM/%.composition.qza: Output/ANCOM/%.qza | build.log
	$(DOCKER) qiime composition add-pseudocount --i-table $(addprefix /,$(word 1,$^)) --o-composition-table $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%.qzv: Output/ANCOM/%.qza Output/Metadata/whole.tsv | build.log
	$(DOCKER) qiime composition ancom --i-table $(addprefix /,$(word 1,$^)) --m-metadata-file $(addprefix /,$(word 2,$^)) --m-metadata-column "Site" --o-visualization $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/ANCOM/%/ancom.tsv: Output/ANCOM/%.qzv | build.log
	$(DOCKER) rm -rfv $(addprefix /,$(@D))
	$(DOCKER) qiime tools export --input-path $(addprefix /,$(word 1,$^)) --output-path $(addprefix /,$(@D))

step08: Output/ANCOM/whole.DADA2.gg.composition/ancom.tsv Output/ANCOM/whole.DADA2.silva.composition/ancom.tsv Output/ANCOM/whole.Deblur.gg.composition/ancom.tsv Output/ANCOM/whole.Deblur.silva.composition/ancom.tsv
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
