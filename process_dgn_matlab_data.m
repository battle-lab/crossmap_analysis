input_geno_mat_fn = '/work-zfs/abattle4/lab_data/dgn/processed_data/eqtl_data/final_completed_genotype.mat'
geno_fn = '/work-zfs/abattle4/ashis/progdata/dgn/genotype/final_completed_genotype.txt'
annot_fn = '/work-zfs/abattle4/ashis/progdata/dgn/genotype/snp_annot.txt'

load(input_geno_mat_fn)

fid = fopen(geno_fn,'wt');

fprintf(fid,'id');
fprintf(fid,'\t%s',completed_genotype.ld_subjects{:});
fprintf(fid,'\n');
for ii = 1:size(completed_genotype.data,2)
  if rem(ii,1000) == 0
    disp(ii)
  end
  fprintf(fid, completed_genotype.snps{ii});
  fprintf(fid,'\t%d',completed_genotype.data(:,ii));
  fprintf(fid,'\n');
end
fclose(fid)


% process snp annotation file
fh = fopen(annot_fn, 'wt');
fprintf(fh, 'snp\tchr\tpos\n');
for(i=1:size(completed_genotype.snps,1))
  if mod(i,10000)==0
display(i);
end
fprintf(fh, '%s\tchr%d\t%d\n', char(completed_genotype.snps(i)), completed_genotype.snp_chrom(i), completed_genotype.snp_pos(i));
end
fclose(fh);
