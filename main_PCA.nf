#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Author: Justin Pelletier (inspired by Peyton McClelland script)
 * Year: 2025
 * Pipeline: LASER PCA Projection with QC 
 */

params.nPCs      = 20
params.min_prob  = 0
params.seed      = 11


process qc_norm_ref {
    tag "$chr"
    input:
        tuple val(chr), path(vcf), path(vcf_tbi)
    output:
        tuple val(chr), path("ref_${chr}.qc.vcf.gz"), path("ref_${chr}.qc.vcf.gz.tbi")
    script:
    """
    module load bcftools
    bcftools norm -m -both -f "${params.ref_fasta}" "$vcf" | \
      bcftools view -f PASS -q 0.05 -Q 0.95 | \
      bcftools annotate -x INFO,^GT | \
      bcftools view -S "${params.qc_ref_list}" --force-samples | \
      bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ref_${chr}.qc.vcf.gz
    tabix -f ref_${chr}.qc.vcf.gz
    """
}

process qc_norm_study {
    tag "$chr"
    input:
        tuple val(chr), path(vcf), path(vcf_tbi)
    output:
        tuple val(chr), path("study_${chr}.qc.vcf.gz"), path("study_${chr}.qc.vcf.gz.tbi")
    script:
    """
    module load bcftools
    bcftools norm -m -both -f "${params.ref_fasta}" "$vcf" | \
      bcftools view -q 0.05 -Q 0.95 | \
      bcftools annotate -x INFO,^GT | \
      bcftools view -S "${params.qc_study_list}" --force-samples | \
      bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o study_${chr}.qc.vcf.gz
    tabix -f study_${chr}.qc.vcf.gz
    """
}

process intersect {
  tag "$chr"
  input:
    tuple val(chr), path(ref), path(ref_tbi), path(study), path(study_tbi)
  output:
    tuple val(chr), path("${chr}.merged.vcf.gz"), path("${chr}.merged.vcf.gz.tbi")
  script:
  """
  module load bcftools
  set -euo pipefail

  tmp=\$(mktemp -d)
  trap 'rm -rf "\$tmp"' EXIT

  bcftools isec -n=2 -p "\$tmp" -Oz "$ref" "$study"
  bcftools merge -m all -Oz -o "${chr}.merged.vcf.gz" "\$tmp/0000.vcf.gz" "\$tmp/0001.vcf.gz"
  tabix -f "${chr}.merged.vcf.gz"
  """
}



process final_qc_and_prune {
  tag "$chr"
  input:
    tuple val(chr), path(vcf), path(vcf_tbi)

  output:
    tuple val(chr), path("${chr}.pruned.vcf.gz"), path("${chr}.pruned.vcf.gz.tbi")

  script:
  """
  module load StdEnv/2020 plink/1.9b_6.21-x86_64 bcftools

  # ─── prepare low complexity bed ────
  # if header_bed==TRUE, drop the first line, then keep only cols 1-3.
  if [ "${params.header_bed}" = "TRUE" ]; then
    tail -n +2 ${params.lowcomplexity_bed} | cut -f1-3 > lowcomplexity.clean.bed
  else
    cut -f1-3 ${params.lowcomplexity_bed} > lowcomplexity.clean.bed
  fi

  # ─── Pre-filter with bcftools to exclude regions ───
  bcftools view -T ^lowcomplexity.clean.bed "$vcf" -Oz -o ${chr}.filtered.vcf.gz

  # ─── Run PLINK on the filtered VCF ───
  plink --vcf ${chr}.filtered.vcf.gz \
        --double-id \
        --real-ref-alleles \
        --snps-only \
        --maf 0.01 \
        --geno 0.001 \
        --biallelic-only strict \
        --make-bed \
        --out ${chr}.final_qc


  # ─── LD pruning ────
  plink --bfile ${chr}.final_qc \
        --indep-pairwise 200 50 0.2 \
        --out ${chr}.prune

  plink --bfile ${chr}.final_qc \
        --extract ${chr}.prune.prune.in \
        --make-bed \
        --out ${chr}.pruned

  # ─── back to VCF ────
  plink --bfile ${chr}.pruned \
        --recode vcf bgz \
        --out ${chr}.pruned

  tabix -f ${chr}.pruned.vcf.gz
  """
}



process concat_pruned_vcfs {
    input: path vcf_files
    output: path("allchr.pruned.vcf.gz")
    script:
    """
    module load bcftools
    bcftools concat -Oz -o allchr.pruned.vcf.gz ${vcf_files.join(' ')}
    """
}

//////////////////////////////
// 1) Double up the sample‐ID lists
//////////////////////////////

process double_ref_ids {
  input:
    path qc_ref_list
  output:
    path "ref_list.id"
  script:
  """
  awk '{print \$1\"_\"\$1}' ${qc_ref_list} > ref_list.id
  """
}

process double_study_ids {
  input:
    path qc_study_list
  output:
    path "study_list.id"
  script:
  """
  awk '{print \$1\"_\"\$1}' ${qc_study_list} > study_list.id
  """
}



//////////////////////////////
// 2) Split the doubled study list into 1 000‐line chunks
//////////////////////////////

process split_study_list {
  input:
    path study_list
  output:
    path "study_batch_*.id"
  script:
  """
  split -d -l ${params.batch_size} --additional-suffix .id ${study_list} study_batch_
  """
}


//////////////////////////////
// 3) Prepare & PCA the *reference* ONCE
//////////////////////////////

process prepare_reference {
  tag "prepare_reference"
  input:
    path vcf
    path ref_list
  output:
    path "reference_pruned.vcf.gz"
    path "reference_pruned.vcf.gz.tbi"
    path "reference_pruned.geno"
    path "reference_pruned.site"
    path "reference.RefPC.coord"
  publishDir "results/LASER_PCA", mode: "copy"
  script:
  """
  module load StdEnv/2020 gcc/9.3.0 bcftools

  # 1) extract only reference samples
  bcftools view -S ${ref_list} --force-samples -Oz -o reference_pruned.vcf.gz ${vcf}
  tabix -f reference_pruned.vcf.gz

  # 2) convert to geno/site
  ${params.path_to_laser}/vcf2geno/vcf2geno --inVcf reference_pruned.vcf.gz --out reference_pruned

  # 3) run PCA on reference
  ${params.path_to_laser}/laser \
    -g reference_pruned.geno \
    -k ${params.nPCs} \
    -pca 1 \
    -o reference
  """
}


//////////////////////////////
// 4) Project each 1000‐sample batch in parallel
//////////////////////////////

process laser_pca_projection_batch {
  publishDir "results/LASER_PCA/batches", mode: "copy"

  input:
   tuple path(batch_file),
          path(vcf),
          path(ref_geno),
          path(ref_site),
          path(ref_pca)
  
  output:
    path "trace_*.ProPC.coord"

  script:
  """
  module load bcftools

  batch_id=\$(basename ${batch_file} .id)

  # extract just this batch
  bcftools view -S ${batch_file} --force-samples -Oz \
    -o study_\${batch_id}.vcf.gz ${vcf}
  tabix -f study_\${batch_id}.vcf.gz

  # convert to geno/site
  ${params.path_to_laser}/vcf2geno/vcf2geno --inVcf study_\${batch_id}.vcf.gz \
                                          --out study_\${batch_id}

  # project onto reference PCs
  cat <<EOF > trace_\${batch_id}.conf
STUDY_FILE    study_\${batch_id}.geno
GENO_FILE     ${ref_geno}
COORD_FILE    ${ref_pca}
FIRST_IND     1
LAST_IND      \$(wc -l < study_\${batch_id}.geno)
OUT_PREFIX    trace_\${batch_id}
DIM           ${params.nPCs}
EOF

  ${params.path_to_laser}/trace -p trace_\${batch_id}.conf
  """
}



//////////////////////////////
// 5) Finally, merge all batches with the reference PCA
//////////////////////////////

process merge_proj_coords {
  input:
    path ref_pca
    path proj_files
  output:
    path "final.ProPC.coord"
  publishDir "results/LASER_PCA", mode: "copy"
  script:
  """
  # ─── 1) Clean reference PCA and write to trimmed output ───
  awk 'NR==1 { print; next }
       {
         # Clean FID and IID: ref1_ref1 → ref1
         for (i=1; i<=2; i++) {
           split(\$i, parts, "_");
           \$i = parts[1];
         }
         # Print all fields with tab separation
         printf "%s\\t%s", \$1, \$2;
         for (i=3; i<=NF; i++) printf "\\t%s", \$i;
         printf "\\n";
       }' ${ref_pca} > final.ProPC.coord

  # ─── 2) Append cleaned batch projections to trimmed output ───
  for f in ${proj_files.join(' ')}; do
    tail -n +2 \$f | awk '{
      # Clean FID and IID: sample1_sample1_sample1_sample1 → sample1_sample1
      for (i=1; i<=2; i++) {
        split(\$i, parts, "_");
        \$i = parts[1] "_" parts[2];
      }
      # Keep only popID, indivID, and PCs (drop L, K, t, Z)
      printf "%s\\t%s", \$1, \$2;
      for (i=7; i<=NF; i++) printf "\\t%s", \$i;
      printf "\\n";
    }' >> final.ProPC.coord
  done

  # ─── 3) Merge batch projections only (full version) ───
  # Write header from first file
  head -n 1 ${proj_files[0]} > final.ProPC.full.coord

  # Append cleaned rows from all batch files
  for f in ${proj_files.join(' ')}; do
    tail -n +2 \$f | awk '{
      # Clean FID and IID: sample1_sample1_sample1_sample1 → sample1_sample1
      for (i=1; i<=2; i++) {
        split(\$i, parts, "_");
        \$i = parts[1] "_" parts[2];
      }
      # Print all fields with tab separation
      printf "%s\\t%s", \$1, \$2;
      for (i=3; i<=NF; i++) printf "\\t%s", \$i;
      printf "\\n";
    }' >> final.ProPC.full.coord
  done
  """
}



process run_pca_analysis {
  tag "run_pca_analysis"

  input:
    path pca_file
  output:
    path "*"
  publishDir "results/PCA_plots", mode: "copy"

  script:
  """
  module load r/4.5.0

  mkdir -p plot_PCA
  Rscript /lustre07/scratch/justinp/Imputation_PAPER/FINAL_2025/PCA/parallel_PCA_nf/bin/PCA.R ${pca_file} ${params.meta_file} ${params.qc_study_list} ${params.k} ${params.n_pcs} ${params.threshold_N}
  """
}




workflow {

    lowcomplexity_bed_ch = Channel.fromPath(params.lowcomplexity_bed)

     ref_ch = Channel.fromPath(params.input_reference, checkIfExists:true)
                    .map { f ->
                        def match = (f.name =~ /chr(\d+)/)
                        if (match) {
                            def chr = match[0][1] as Integer
                            if (chr >= 1 && chr <= 22)
                                tuple(chr.toString(), f, file(f.toString() + '.tbi'))
                        }
                    }
                    .filter { it != null }

    study_ch = Channel.fromPath(params.input_study, checkIfExists:true)
                      .map { f ->
                          def match = (f.name =~ /chr(\d+)/)
                          if (match) {
                              def chr = match[0][1] as Integer
                              if (chr >= 1 && chr <= 22)
                                  tuple(chr.toString(), f, file(f.toString() + '.tbi'))
                          }
                      }
                      .filter { it != null }


    //ref_ch.view { "REF → $it" }
    //study_ch.view { "STUDY → $it" }

    // turn those two params into channels
    qc_ref_list_ch   = Channel.fromPath(params.qc_ref_list)
    qc_study_list_ch = Channel.fromPath(params.qc_study_list)

    // Normalization, QC and filtering
    ref_qc   = qc_norm_ref(ref_ch)
    study_qc = qc_norm_study(study_ch)

    // Join outputs by chr
    ref_study_joined = ref_qc.join(study_qc)

    // Intersect shared variants
    intersect_results = intersect(ref_study_joined)

    // Final QC + LD pruning
    pruned_all = final_qc_and_prune(intersect_results)

    //
    // 1) concatenate all the per-chr pruned VCFs
    //
    concat_results = concat_pruned_vcfs(
      pruned_all.map { it[1] }.collect()
    )

    //concat_results.view { "CONCAT VCF → $it" }

    //
    // 2) build the doubled ID lists
    //
    ref_ids   = double_ref_ids(qc_ref_list_ch)
    study_ids = double_study_ids(qc_study_list_ch)


    //
    // 3) prepare the reference PCA once (emit: ref_vcf, ref_tbi, ref_geno, ref_site, ref_pca)
    //
    ( ref_vcf, ref_tbi, ref_geno, ref_site, ref_pca ) =  prepare_reference(concat_results, ref_ids)

    ref_geno.view { "REF GENO → $it" }
    ref_site.view { "REF SITE → $it" }
    ref_pca.view { "REF PCA → $it" }


    //
    // 4) split the study IDs into 1 000-sample chunks
    //
    study_batches = split_study_list(study_ids)
      .flatMap { it }    
    study_batches.view { "STUDY_BATCHES → $it" }


    // 5) project each batch in parallel, feeding all five channels
    proj_input = study_batches
      .combine(concat_results)
      .combine(ref_geno)
      .combine(ref_site)
      .combine(ref_pca)
      .map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) }


    proj_input.view { "PROJ_INPUT → ${it[0].name}" }

    projected = laser_pca_projection_batch(proj_input)

    projected.collect().view { "MERGE INPUT FILES → $it" }


    // 7) finally stitch them back together
    //final_pca = merge_proj_coords(ref_pca, projected)
    final_pca = merge_proj_coords(ref_pca, projected.collect())



    // 8 Generates the plots
    pca_plots = run_pca_analysis(final_pca)
}

