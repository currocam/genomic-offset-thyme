# This script creates an Rds file from a VCF file and a text file from SLiM
# It doesn't follow any standard format (sadly I didn't find any)

using VCFTools, RCall
infilevcf = popfirst!(ARGS)
infiletxt = popfirst!(ARGS)
outfile = popfirst!(ARGS)

function main(infilevcf::String, infiletxt::String, outfile::String)
  # Read metadat from the txt file
  simulation = Dict()
  file = open(infiletxt, "r")
  while true
          keys = readline(file)
          if isempty(keys)
                  break
          end
          simulation[keys[2:end]] = parse.(Float64, split(readline(file)))
  end
  close(file)
  # Read genotypes from the VCF file
  simulation["Genotype"] = convert_gt(
          Int8, infilevcf; model=:additive,
          impute=false, center=false, scale=false
  )
  # Match the SNPs from the VCF file with the SNPs from the txt file
  raw_pos = read(`bash -c "grep -v '#' < $infilevcf | cut -f 2"`, String)
  pos = parse.(Int, split(raw_pos))
  position2index =  x -> [findfirst(y -> y == snp, pos) for snp in x]
  trait_index = 1
  while true
          trait = "QTLs trait $trait_index"
          if !haskey(simulation, trait)
                  break
          end
        snps = Int.(simulation[trait])
        indexes = position2index(snps)
        simulation[trait] = Int.(snps[indexes .!== nothing])
        simulation["Index QTLs $trait_index"] = Int.(indexes[indexes .!== nothing])
        trait_index += 1
  end
  # Convert simulation dictionary to R object
  @rput simulation
  # Save as compress Rds file
  R"saveRDS(simulation, $outfile, compress = 'xz')"
end

main(infilevcf, infiletxt, outfile)

