#nohup perl ../DMG_upstream2k_gene_downstram2k.pl ./filter_depth4_nC5.dmr ichang_chromosome_filter267+234gene.sort.gff3 filter_depth4_nC5_gene_CG.dmr filter_depth4_nC5_upstream_2k_CG.dmr filter_depth4_nC5_downstream_2k_CG.dmr&
($in1,$in2,$out1,$out2,$out3)=@ARGV;
open IN1,"$in1";
open IN2,"$in2";
open OUT1,">$out1";
open OUT2,">$out2";
open OUT3,">$out3";
while(chomp($line=<IN1>)){
	@tmp=split(" ",$line);
	$chr=$tmp[0];
    $start=$tmp[1];
	$id=$chr."_".$start;
       	$hash{$id}=$line;;
}
while (chomp($line=<IN2>)) {
    @tmp=split(" ",$line);
    $chr=$tmp[0];
    $feature=$tmp[2];
    $strand=$tmp[6];
    $start=$tmp[3];
    $end=$tmp[4];
    if ($feature eq "gene") {
        $mark=0;
        if ($line=~/(.*)ID=(.*?);(.*)/) {
            $gene=$2;
						push(@gene,$gene);
        }
				foreach($start..$end){
					$id=$chr."_".$_;
					if(exists($hash{$id})){
						print OUT1 "$gene\t$start\t$end\t$hash{$id}\n";
					}
        }
            if ($strand eq "+") {
                $start_promoter=$start-2001;
                $end_promoter=$start-1;
            }
			if ($strand eq "-") {
                $start_promoter=$end+1;
                $end_promoter=$end+2001;
            }
						foreach($start_promoter..$end_promoter){
							$id=$chr."_".$_;
							if(exists($hash{$id}))
							{
								print OUT2 "$gene\t$start_promoter\t$end_promoter\t$hash{$id}\n";
							}
						}
            if ($strand eq "-") {
                $start_downstream=$start-2001;
                $end_downstream=$start-1;
            }
			if ($strand eq "+") {
                $start_downstream=$end+1;
                $end_downstream=$end+2001;
            }
						foreach($start_downstream..$end_downstream){
							$id=$chr."_".$_;
							if(exists($hash{$id}))
							{
								print OUT3 "$gene\t$start_downstream\t$end_downstream\t$hash{$id}\n";
							}
						}	
						
    }
}