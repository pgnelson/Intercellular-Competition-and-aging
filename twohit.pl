#!/usr/bin/perl
use strict;
use Math::Random qw(random_normal random_multinomial random_exponential random_poisson);
#use Math::BigFloat;
#	Math::BigFloat->precision(-1);
#	Math::BigFloat->accuracy(2);

#use Bignum;
my $num_vloci = 1;
my $num_yloci = 1;
my $num_cells = 10**5;
my $prob_mut = 10**(-6);
my $prob_v_mut = 10**-1;#$num_vloci * $prob_mut * $num_cells;
my $prob_y_mut = 10**-1;#$num_yloci * $prob_mut * $num_cells;
my $Ne = 10000;
my $muv = -1*10**-2;#$prob_v_mut / ($num_vloci * $num_cells);
my $muy = -1*10**-2;#$prob_y_mut / ($num_yloci * $num_cells)*10**-3;
my $c = 0.5;
my $r_max;
my $generations = 10000;
my $max_lifespan = 10000;
my $print_within_life = 0;
my $evolvingmu = 0;
my $within_lifespan_num_printouts = 1000;
my $malthusian = 1;
my $birth_scalar = 0.001;
my $frac_z_affects_death = 1;
#my @cancer_stages = (1, 0, 0);
#my @vigor_stages = (1, 0, 0);
open (OUTPUT, '>comp.txt');
#printf OUTPUT "replicate\tlog(muv)\tlog(muy)\tlog(r)\tfitness\n";
my $death_cutoff = 0.0001;

my $optimizing_r = 1;
my $get_life_history =0;

printf OUTPUT "survival\tlifetime_fitness\tscaler\tr_max\tmuv\tmuy\text_death\tvig_2nd_stage\tcan_2nd_stage\trobustness\ttime\tabs_time\torg_fit\t";
printf OUTPUT "m1_1	m1_2	m1_3	m2_1	m2_2	m2_3	m3_1	m3_2	m3_3 \n";


for (my $parameter1 = 0; $parameter1 <= 1; $parameter1+=1){	
	my $start = 0.01;#0.00183645625654775 - 0.0001;
	my $fin = 0.8;#0.00183645625654775 + 0.0001;
	for (my $parameter2 = 0; $parameter2 <= 1; $parameter2+=1){
		my $num_cancer_death = 0;
		my $ext_death = 0;#$parameter1;
		my $muv = -10**(-3);
		my $muy = -10**(-5);
		my $r = $parameter1;#0.00673288099381955;

		my @vigor_stages = (1, 0, 0);
		my @cancer_stages = (1, 0, 0);

		if ($optimizing_r == 1){
			print "death_cutoff\tfrac_z_death\text_death\tmuv\tmuc";
			printf OUTPUT "death_cutoff\tfrac_z_death\text_death\tmuv\tmuc\tv2ndstg\tc2ndstg\ttime\tlifetime_fitness\tmean_vigor\tmean_coop\tscaler\tcell_comp";
			for (my $i = 0; $i < scalar(@vigor_stages); $i++){
				for (my $j = 0; $j < scalar(@cancer_stages); $j++){
					print "\t" . $vigor_stages[$i] . "_" . $cancer_stages[$j];
					printf OUTPUT "\t" . $vigor_stages[$i] . "_" . $cancer_stages[$j];
				}
			}
			print "\n";
			printf OUTPUT "\n";
		}

		for (my $replicate = 0; $replicate <= 0.0; $replicate+=0.0001){
			$death_cutoff = 0.5;#$replicate;
			my $vigor_2nd_stage = $parameter1;
			my $cancer_2nd_stage = $parameter2;
			$r=$replicate;
			#$num_cancer_death = 0;		
			$malthusian = 1; #$replicate;
			$frac_z_affects_death = 0;#$parameter1;
			if	($get_life_history == 1){
				$print_within_life = 0;
				my $scaler = find_scaler($r, $muv, $muy, 0.5, $max_lifespan, $ext_death, $frac_z_affects_death, $vigor_2nd_stage, $cancer_2nd_stage);	
				$scaler = find_scaler($r, $muv, $muy, $scaler, $max_lifespan, $ext_death, $frac_z_affects_death, $vigor_2nd_stage, $cancer_2nd_stage);	

				my @fitness = @{simulate_lifespan($r, $muv, $muy, $scaler, $ext_death, $print_within_life, $malthusian, $frac_z_affects_death, $vigor_2nd_stage, $cancer_2nd_stage)};
				print $death_cutoff . "\t" .$frac_z_affects_death . "\t". $ext_death . "\t" . log(-$muv)/log(10) . "\t" .log(-$muy)/log(10) . "\t" . $vigor_2nd_stage. "\t" .  $cancer_2nd_stage;
				printf OUTPUT $death_cutoff . "\t" .$frac_z_affects_death . "\t". $ext_death . "\t" . $r . "\t" . log(-$muv)/log(10) . "\t" .log(-$muy)/log(10) . "\t" . $vigor_2nd_stage. "\t" .  $cancer_2nd_stage;
				foreach(@fitness){
					print "\t" . $_;
					printf OUTPUT "\t" . $_;
				}
				print "\n";
				printf OUTPUT "\n";
			}
			if ($optimizing_r == 1){
				my @optimal = @{optimize_r($muv, $muy, $ext_death, $frac_z_affects_death, $vigor_2nd_stage, $cancer_2nd_stage)};
				print $death_cutoff . "\t" .$frac_z_affects_death . "\t". $ext_death . "\t" . log(-$muv)/log(10) . "\t" .log(-$muy)/log(10) . "\t" . $vigor_2nd_stage. "\t" .  $cancer_2nd_stage;
				printf OUTPUT $death_cutoff . "\t" .$frac_z_affects_death . "\t". $ext_death . "\t" . log(-$muv)/log(10) . "\t" .log(-$muy)/log(10) . "\t" . $vigor_2nd_stage. "\t" .  $cancer_2nd_stage;
				
				foreach(@optimal){
					print "\t" . $_;
					printf OUTPUT "\t" . $_;
				}
				print "\n";
				printf OUTPUT "\n";
			}
		}
		#print $r . "\t" . $num_cancer_death . "\n";
	}
}
close OUTPUT;


sub simulate_lifespan{
	my $r_max = $_[0];
	my $org_muv = $_[1];
	my $org_muy = $_[2];
	my $scaler= $_[3];
	my $ext_death = $_[4];
	my $printout = $_[5];
	my $need_malthusian = $_[6];
	my $frac_death = $_[7];
	my $vig_2nd_stage = $_[8];
	my $can_2nd_stage = $_[9];
	
	
	my @vigor_stages = (1, $vig_2nd_stage, 0);
	my @cancer_stages = (1, $can_2nd_stage, 0);

	my $num_cancer_stages = scalar(@cancer_stages);
	my $num_vigor_stages = scalar(@vigor_stages);
#	print $r_max . "r\t" . $org_muv. "v\t" . $org_muy. "y\t" . $scaler ."s\t";
	my $mean_vigor = 0;
	my $mean_coop = 0;
	my $last  = 0; 
	my @vloci;
	my @yloci;
	my $lifetime_fitness = 0;
	my @cell_type = []; 
	for (my $i = 0; $i < $num_vigor_stages; $i++){
		for (my $j = 0; $j < $num_cancer_stages; $j++){
			$cell_type[$i][$j] = 0;
		}
	}
	$cell_type[0][0] = 1;

	my $time = 0;
	my $vigor = 0;
	my $cancer = 0;
	
	my $surviving = 1;
	my @inst_fit_list = [];
	$inst_fit_list[0] = 0;
	my $survival = 1;

	while ($last ==0){		
		my $mean_cell_comp = 0;#calculate mean fitness
		for (my $i = 0; $i <= $num_vigor_stages; $i++){
			for (my $j = 0; $j <= $num_cancer_stages; $j++){
				$mean_cell_comp += $cell_type[$i][$j] * (1 - $c * $cancer_stages[$j] ) * $vigor_stages[$i];
				#print $i . "\t" . $j . "\t" . $num_cancer_stages . "\t".$num_vigor_stages ."\t". $mean_cell_comp  . "\t" . $cancer_stages[$j] . "\t" . $vigor_stages[$i] ;
				#<STDIN>;
			}
		}		
		#increment cell classes
		my $dist_mean = 10000;
		for (my $i = 0; $i < $num_vigor_stages; $i++){
			for (my $j = 0; $j < $num_cancer_stages; $j++){
				$cell_type[$i][$j] += $scaler * ($r_max * $cell_type[$i][$j] * (1 - $c * $cancer_stages[$j]) * $vigor_stages[$i] / $mean_cell_comp + (1 - $r_max)*$cell_type[$i][$j] - $cell_type[$i][$j]);
			}
		}


		for (my $i = $num_vigor_stages - 1; $i >= 0; $i--){
			for (my $j = $num_cancer_stages - 1; $j >= 0; $j--){
				#my $curr_muv= -(1/ $dist_mean * $cell_type[$i - 1][$j]) * random_poisson(1, -1*$org_muv * $dist_mean * $cell_type[$i - 1][$j]);
				#my $curr_muy = -(1/ $dist_mean * $cell_type[$i][$j - 1]) * random_poisson(1, -1*$org_muy * $dist_mean * $cell_type[$i][$j - 1]);
				my $curr_mu_vig = 0;
				my $curr_mu_coop = 0;
				if ($i > 0){
					$curr_mu_vig = $scaler * $org_muv * $cell_type[$i - 1][$j];
					$cell_type[$i - 1][$j] +=  $curr_mu_vig;
				}
				if ($j > 0){
					$curr_mu_coop = $scaler * $org_muy * $cell_type[$i][$j - 1];
					$cell_type[$i][$j - 1] += $curr_mu_coop;
				}
				$cell_type[$i][$j] -= $curr_mu_vig + $curr_mu_coop;
			}
		}

		$mean_vigor = 0;
		$mean_coop = 0;
		my $org_fit = 0;
		for (my $i = 0; $i < $num_vigor_stages; $i++){
			for (my $j = 0; $j < $num_cancer_stages; $j++){
				$mean_vigor += $cell_type[$i][$j] * $vigor_stages[$i];
				$mean_coop += $cell_type[$i][$j] * $cancer_stages[$j];
				$org_fit += $cell_type[$i][$j] * $vigor_stages[$i] * $cancer_stages[$j];
			}
		}

	#	if ($mean_vigor >= $death_cutoff && $mean_coop >= $death_cutoff){
		if ($org_fit >= $death_cutoff ){
			$survival *= ((1-$ext_death)*((1-$frac_death)*1 + $frac_death * $org_fit))**($scaler * $birth_scalar);
		}else{
			$survival = 0;
		}
		$lifetime_fitness +=  $survival * (1 * $frac_death + (1 - $frac_death) * $org_fit)* $scaler * $birth_scalar;
		if ($survival <=0 || $time == $max_lifespan - 1){
			$last  = 1; 	
		}elsif($need_malthusian == 1){
			$inst_fit_list[$time] = $org_fit;
		}

		if ($printout == 1 && ($time+1) % ($max_lifespan / $within_lifespan_num_printouts) == 0 ){
			print $survival . "\t" . $lifetime_fitness ."\t" . $scaler . "\t". $r_max . "\t" . $org_muv ."\t" . $org_muy ."\t" . $ext_death . "\t" . $time . "\t" . $org_fit;
			printf OUTPUT $survival . "\t" . $lifetime_fitness ."\t" . $scaler . "\t". $r_max . "\t" . $org_muv ."\t" . $org_muy ."\t" . $ext_death . "\t" .$vig_2nd_stage . "\t" . $can_2nd_stage . "\tm" . $vig_2nd_stage . "_" . $can_2nd_stage. "\t" . $time . "\t" . $time * $scaler ."\t" . $org_fit;
			for (my $i = 0; $i < $num_vigor_stages; $i++){
				for (my $j = 0; $j < $num_cancer_stages; $j++){
					print "\t" . $cell_type[$i][$j];
					printf OUTPUT "\t" . $cell_type[$i][$j];
					$org_fit += $cell_type[$i][$j] * $vigor_stages[$i] * $cell_type[$i][$j] * $cancer_stages[$j];
				}
			}
			print "\n";
			print OUTPUT "\n";
		}
		$time++;
	}
	my $report_scaler = $scaler;
	if ($need_malthusian == 1){
		my $malth_r = 10**-10;
		my $adj_fit_total = 0;
		$time = 0;
		
		$survival = 1;
		foreach(@inst_fit_list){
			if ($_ > 0){
				$survival *= ((1-$ext_death)*((1-$frac_death)*1 + $frac_death * $_))**( $scaler * $birth_scalar );
			}else{
				$survival = 0;
			}
			$adj_fit_total +=  exp(-$time * $malth_r)*$survival * (1 * $frac_death + (1- $frac_death) * $_ )* $scaler * $birth_scalar;
			$time += $scaler * $birth_scalar;
			#((1-exp(-$scaler))/($malth_r*$scaler)) * <-adjustment to for mean exp term
		}		
		my $target;
		my $step;
		my $sign;
		if ($adj_fit_total < 1){
			$target = -1;
			$step = -10**-4;
			$sign = -1;
		}else{
			$target = 1;
			$step = 10**-4;
			$sign = 1
		}
		my $lower_search_bounds = -2;
		my $upper_search_bounds = 2;
		my $min = $lower_search_bounds;
		my $max = $upper_search_bounds;
		#while($sign * $adj_fit_total > $target){
		while(abs($adj_fit_total -1) > 10**-10){
			$adj_fit_total = 0;
			$time = 0;
			$malth_r = ($min + $max)/2;
			$survival = 1;
			foreach(@inst_fit_list){
				if ($_ > 0){
					$survival *= ((1-$ext_death)*((1-$frac_death)*1 + $frac_death * $_))**($scaler * $birth_scalar);
			#print $survival . "\t";

				}else{
					$survival = 0;
				}
				$adj_fit_total += exp(-$time * $malth_r)*$survival * (1 * $frac_death + (1- $frac_death) * $_)* $scaler * $birth_scalar;
				$time += $scaler * $birth_scalar;
			}
			if ($adj_fit_total > 1){
				$min = $malth_r;
			}else{
				$max = $malth_r;
			}
		#	printf OUTPUT $scaler . "\t" . $r_naught. "\t" . $adj_fit_total. "\t" . $time . "\n";
		#	print "rsearch" . "\t" . $scaler . "\t" . $org_muy . "\t". $org_muv . "\t" . $min . "\t" . $max . "\t" .$malth_r ."\t" . $adj_fit_total . "\n";
			#$malth_r += $step;
			if ($malth_r <= $lower_search_bounds){
				$lower_search_bounds-=10;
				$min=$lower_search_bounds;

			}	
			if ($malth_r >= $upper_search_bounds){
				$upper_search_bounds+=10;
				$max+=$upper_search_bounds;
			}			
		}
		$lifetime_fitness = $malth_r;
	}
	#print $time . "\n";
	my @return_array = ($time, $lifetime_fitness, $mean_vigor, $mean_coop, $scaler, $r_max);
	for (my $i = 0; $i < $num_vigor_stages; $i++){
		for (my $j = 0; $j < $num_cancer_stages; $j++){
			push @return_array, $cell_type[$i][$j];
		}
	}

	return \@return_array;
}

sub find_scaler {
	my $r_max = $_[0];
	my $org_muv = $_[1];
	my $org_muy = $_[2];
	my $best_scaler = $_[3];
	my $max_lifespan = $_[4];
	my $ext_death = $_[5];
	my $frac_death = $_[6];
	my $vig_2nd_stage = $_[7];
	my $can_2nd_stage = $_[8];
	
	my @lifespan_result = @{simulate_lifespan($r_max, $org_muv, $org_muy , $best_scaler, $ext_death, 0, 0, $frac_death, $vig_2nd_stage, $can_2nd_stage)};
	my $best_lifespan = $lifespan_result[0];
	#print $best_scaler. "\t" . $best_lifespan . "\t" . $max_lifespan . "\t";
	while($best_lifespan >= $max_lifespan){
		$best_scaler *=1.1;
#				print $best_scaler . "\tscal\n";

		@lifespan_result = @{simulate_lifespan($r_max, $org_muv, $org_muy , $best_scaler, $ext_death, 0, 0, $frac_death, $vig_2nd_stage, $can_2nd_stage)};

		$best_lifespan= $lifespan_result[0];
	#print $best_scaler. "\t" . $best_lifespan . "\tscal\n";
	}
	$best_scaler = $best_scaler *  $best_lifespan / (0.3*$max_lifespan);
	#	print "scalfin\n";
	return $best_scaler;
}

sub optimize_r {
	my $last = 0;
	my $r_curr = 0;
	my $best_r = 0;
	my $muv_curr = $_[0];
	my $muy_curr = $_[1];
	my $ext_death_curr = $_[2];
	my $frac_death = $_[3];
	#print $frac_death . "\t" . $muv_curr . "\t" . $muy_curr . "\t" . $r_curr . "\n";
	my $vig_2nd_stage = $_[4];
	my $can_2nd_stage = $_[5];

	my $a = 0;
	my $b = 1;
	my $phi = (1 + sqrt(5))/ 2;
	my $curr_scaler = find_scaler($a, $muv_curr, $muy_curr, 0.5, $max_lifespan, $ext_death_curr, $frac_death, $vig_2nd_stage, $can_2nd_stage);	

	while (abs($a - $b) > 10**-8){
		my $c = $b + ($a - $b)/$phi;
		my $d = $a + ($b - $a)/$phi;
		$curr_scaler = find_scaler($c, $muv_curr, $muy_curr, $curr_scaler, $max_lifespan, $ext_death_curr, $frac_death, $vig_2nd_stage, $can_2nd_stage);	
		my $fitness_c = ${simulate_lifespan($c, $muv_curr, $muy_curr, $curr_scaler, $ext_death_curr, $print_within_life, $malthusian, $frac_death, $vig_2nd_stage, $can_2nd_stage)}[1];
		$curr_scaler = find_scaler($d, $muv_curr, $muy_curr, $curr_scaler, $max_lifespan, $ext_death_curr, $frac_death, $vig_2nd_stage, $can_2nd_stage);	
		my $fitness_d = ${simulate_lifespan($d, $muv_curr, $muy_curr, $curr_scaler, $ext_death_curr, $print_within_life, $malthusian, $frac_death, $vig_2nd_stage, $can_2nd_stage)}[1];
		if($fitness_c > $fitness_d){
			$b = $d;
		}else{
			$a = $c;
		}
		#print $a . "\t" . $b . "\t" . $c . "\t" . $d . "\t" . $fitness_c . "\t" . $fitness_d . "\n";
	}
	$curr_scaler = find_scaler(($a+$b)/2, $muv_curr, $muy_curr, $curr_scaler, $max_lifespan, $ext_death_curr, $frac_death, $vig_2nd_stage, $can_2nd_stage);
	my @return_array = @{simulate_lifespan(($a+$b)/2, $muv_curr, $muy_curr, $curr_scaler, $ext_death_curr, $print_within_life, $malthusian, $frac_death, $vig_2nd_stage, $can_2nd_stage)};
	#print $return_array[0] . "\t" . $return_array[1] . "\t" . $return_array[2] ."\n";
	return \@return_array;
}
close OUTPUT;
