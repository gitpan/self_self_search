#!/usr/bin/perl
# Last Update by /gn0/jong/Perl/update_subroutines.pl: Fri Sep 19 03:30:56 BST 1997
#_____________________________________________________________________
# Title     : self_self_search.pl
# Usage     : self_self_search.pl <FASTA_db.fa> [o c s m r]
# Function  : self_to_self input database search. This matches all the
#              sequences to others in the given database(sequence)
#              Default search program is 'fasta'. You can specify it
#              to ssearch by using 'a=ssearch' prompt option.
#             It will produce 2 char size subdirectories. In them, you will
#              get search results stored in MSP (matching sequence pair)
#              file format.
# Example   : self_self_search.pl PDB_db.fa
# Keywords  : self_to_self_db_search, self self search, database search
# Options   :
#             DB=   for target DB  "DB=MY_FASTA.fa"
#             File= to get file base(root) name.  "File=TARGET_DB.fa"
#             m  for MSP format directly from FASTA or Ssearch result than through sso_to_msp to save mem
#             s  for the big single output (msp file output I mean)
#             o  for overwrite existing xxxx.fa files for search
#             c  for create SSO file (sequence search out file)
#             r  for reversing the protein sequence !!
#             R  for attaching ranges of sequences
#             k= for k-tuple value. default is 1 (ori. FASTA prog. default is 2)
#             u= for $upper_expect_limit  ## The higher the slower
#             l= for $lower_expect_limit
#             a= for choosing either fasta or ssearch algorithm
#             d= for defining the size of subdir made. 2 means it creates
#                    eg, DE while 1 makes D
#             d  for $make_gz_in_sub_dir_opt, putting resultant sso files in gz format and in single char subdir
#             D  for $make_msp_in_sub_dir_opt, convert sso to msp and put in sub dir like /D/, /S/
#             n  for new format (msp2 format)
#
#  $over_write     = o by o -o
#  $create_sso     = c by c -c
#  $single_big_msp = s by s -s
#  $msp_directly_opt = m by m -m
#  $reverse_query    =r by r -r
#  $upper_expect_limit = by u=
#  $machine_readable = M by M -M
#  $do_in_batch = b by b -b
#  $make_msp_in_sub_dir_opt = D by D -D
#  $make_gz_in_sub_dir_opt = d by d -d
#  $new_format=n by n -n
#  $sub_dir_size= by d=
#  $add_range=R by R -R
#
# Returns   :
# Argument  :
# Version   : 1.3
#-----------------------------------------------------------------------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# All the defaults
#______________________________
$sub_dir_size=2;
$make_msp_in_sub_dir_opt='D';
$machine_readable='M';

$pdb_seq_fasta    =$ENV{'PDB_SEQ_FASTA'};  ## default pdb is p100
$pdb40_seq_fasta  =$ENV{'PDB40_SEQ_FASTA'};  ## default pdb is p100
$pdb100_seq_fasta =$ENV{'PDB100_SEQ_FASTA'};  ## default pdb is p100
$pdbd100_seq_fasta=$ENV{'PDBD100_SEQ_FASTA'};  ## default pdb is p100
$pdbd40_seq_fasta =$ENV{'PDBD40_SEQ_FASTA'};  ## default pdb is p100
$owl_db_fasta     =$ENV{'OWL_FASTA'};
$pdb_seq          =$ENV{'PDB_SEQ'};
$pdbd_seq         =$ENV{'PDBD_SEQ'};
$swiss_db_fasta   =$ENV{'SWISS_FASTA'};
$pdbd95_seq_fasta =$ENV{'PDB95D_FASTA'};
$pdb40d_old_fasta =$ENV{'PDB40D_OLD_FASTA'};  ## default pdb is p100
$all_bac_genome   =$ENV{'ALL_BAC_GENOME'};
$pdb95d_interm_lib=$ENV{'PDB95D_INTERM_LIB'};


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
# preprocessing the inputs
#______________________________
@file=@{&parse_arguments(1)};
@file=@{&scramble_array(\@file)};

&show_options;

unless($upper_expect_limit=~/\S/){
   $upper_expect_limit=10;
}

&self_self_search(\@file, $over_write, $msp_directly_opt, "u=$upper_expect_limit", $do_in_batch,
                   $create_sso, $reverse_query, $single_big_msp, $machine_readable,
                   "DB=$owl_db_fasta", $make_msp_in_sub_dir_opt, "d=$sub_dir_size",
				   $new_format,   $add_range);

#______________________________________________________________________________
# Title     : self_self_search
# Usage     : &self_self_search(\@file, $over_write, $msp_directly_opt, $create_sso, $single_big_msp);
# Function  : self_to_self input database search with reverse query as an option
# Example   : &self_self_search(\@file, $over_write, $msp_directly_opt, $create_sso, $single_big_msp);
# Warning   :
# Keywords  : do_self_self_search, self_self_sequence_search, self_self_seq_search,
#             self_to_self_search, self_to_rev_self_search, self_to_reversed_self_search
# Options   :
#             Query_seqs=  for enquiry sequences eg)  "Query_seqs=$ref_of_hash"
#             DB=   for target DB  "DB=$DB_used"
#             File= to get file base(root) name.  "File=$file[0]"
#             m  for MSP format directly from FASTA or Ssearch result than through sso_to_msp to save mem
#             s  for the big single output (msp file output I mean)
#             o  for overwrite existing xxxx.fa files for search
#             c  for create SSO file (sequence search out file)
#             r  for reverse the query sequence
#             R  for attaching ranges of sequences
#             b  for doing in batch. Reads all the seqs in memory at one time
#             m10 for machine readable form
#             k= for k-tuple value. default is 1 (ori. FASTA prog. default is 2)
#             u= for $upper_expect_limit
#             l= for $lower_expect_limit
#             a= for choosing either fasta or ssearch algorithm
#             d= for defining the size of subdir made. 2 means it creates
#                    eg, DE while 1 makes D
#             d  for $make_gz_in_sub_dir_opt, putting resultant sso files in gz format and in single char subdir
#             D  for $make_msp_in_sub_dir_opt, convert sso to msp and put in sub dir like /D/, /S/
#             n  for new format (msp2 format)
# Version   : 1.9
#-------------------------------------------------------------------------------
sub self_self_search{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	my (%fasta_seqs, $new_format, $over_write, $upper_expect_limit, $lower_expect_limit,
	    $Score_thresh, $margin, $single_big_msp, $sequence_DB, $create_sso, $k_tuple,
        $msp_directly_opt, $machine_readable, $make_gz_in_sub_dir_opt, $sub_dir_size,
	    $do_in_batch, $make_msp_in_sub_dir_opt, $new_format, $machine_readable,
	    $add_range, $reverse_sequence, $sub_dir_head );
	my $algorithm='/gn0/sat/Apps/fasta3_non_shared/fasta3';
	my $msp_directly_opt='m';
    $sub_dir_size=2;          # default
	#$single_big_msp ='s';
    #$create_sso='c';
	$upper_expect_limit=2;
    $k_tuple=1;

	if($vars{'a'}=~/\S+/){ $algorithm          = $vars{'a'}            };
	if($vars{'u'}=~/\d+/){ $upper_expect_limit = $vars{'u'}            };
	if($vars{'l'}=~/\d+/){ $lower_expect_limit = $vars{'l'}            };
	if($vars{'k'}=~/\d+/){ $k_tuple            = $vars{'k'}            };
	if($vars{'t'}=~/\d+/){ $Score_thresh       = $vars{'t'}            };
    if($vars{'m'}=~/\d+/){ $margin             = $vars{'m'}            };
    if($vars{'d'}=~/\d+/){ $sub_dir_size       = $vars{'d'}            };
	if($vars{'s'}=~/\S+/){ $single_big_msp     = 's'                   };
	if($vars{'DB'}=~/\S+/){ $sequence_DB       = $vars{'DB'}           };
	if($vars{'File'}=~/\S+/){ $input_file_name = $vars{'File'}         };
	if($vars{'Query_seqs'}=~/\S+/){ %seq_input = %{$vars{'Query_seqs'}}};
	if($vars{'u'}         =~/\S+/){ $E_val     = $vars{'u'}            };
	if($char_opt=~/R/){    $add_range          = 'r' }
	if($char_opt=~/o/){    $over_write         = 'o' }
    if($char_opt=~/c/){    $create_sso         = 'c' }
	if($char_opt=~/s/){    $single_big_msp     = 's'; print "\n# Single file opt is set\n"; }
	if($char_opt=~/m/){    $msp_directly_opt   = 'm' }
	if($char_opt=~/M/){    $machine_readable   = 'M' }
	if($char_opt=~/d/){$make_gz_in_sub_dir_opt = 'd' } # for simple search and storing in gz file (sso file will be zipped
	if($char_opt=~/D/){$make_msp_in_sub_dir_opt= 'D' } # for simple search and storing msp file
 	if($char_opt=~/b/){    $do_in_batch        = 'b' } # for reading in all the
    if($char_opt=~/n/){    $new_format         = 'n' }
    if($char_opt=~/r/){ $reverse_sequence   = 'r'  };

    if($do_in_batch=~/b/){
	   for($i=0; $i< @file; $i++){
		  my $input_db_file=$file[$i];
		  %fasta_seqs=%{&open_fasta_files(\$input_db_file)};
		  if($reverse_sequence){ ## reverse the query seqs.
			 %fasta_seqs=%{&reverse_sequences(\%fasta_seqs)};
		  }
          print "\n# self_self_search : \$do_in_batch is set with DB=$input_db_file, File=$input_db_file\n";
		  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		  #  Main sequence search
		  #___________________________________________________
		  my @file_created=@{&do_sequence_search(\%fasta_seqs,
                                                 "d=$sub_dir_size",
                                                 "DB=$input_db_file",
                                                 "File=$input_db_file",
                                                 $single_big_msp,
                                                 $over_write,
                                                 "u=$upper_expect_limit",
                                                 "l=$lower_expect_limit",
                                                 "k=$k_tuple",
                                                 $add_range,
                                                 $create_sso,
                                                 "t=$Score_thresh",
                                                 "m=$margin",
                                                 "a=$algorithm",
                                                 $msp_directly_opt,
                                                 $machine_readable,
                                                 $new_format )};
		  print "\n# File created: @file_created \n";
	   }
   }else{ ## reads in the big database file continously
	   my $make_gz_in_sub_dir_opt='d';
	   my ($ori_seq_name, $first_char, $seq_file_msp_name,  $seq, $seq_name, $first_char);
	   for($i=0; $i< @file; $i++){
	       unless(-s  $file[$i]){
              if(-s "../$file[$i]"){    $input_db_file="../$file[$i]";
              }elsif( -s "../../$file[$i]"){  $input_db_file="../../$file[$i]";
                   print "\n# self_self_search : I found $file[$i] at ../../ ";
              }
		   }else{
              print "\n# self_self_search : I found $file[$i] " if $debug;
		   }

		   open(FASTA, "$file[$i]");
		   while(<FASTA>){
			  if( /\> *((\w\S)\S*)/ ){
				  $ori_seq_name=$1;
				  if($seq=~/\S/ and $seq_name=~/\S/){
                     $seq_file_name="$seq_name\.fa";
                     $seq_file_msp_name="$seq_name\.msp";
                     $seq_file_msp_gz_name="$seq_name\.msp\.gz";
                     $first_char=substr("\U$seq_name", 0, $sub_dir_size);
                     if( $over_write !~/o/ and (-s "$first_char\/$seq_file_msp_name" or -s "$first_char\/$seq_file_msp_gz_name") ){
						 print "\n# $first_char\/$seq_file_msp_name already exists ";
						 $seq='';
					 }else{
                         if($reverse_sequence){ ## reverse the query seqs.
                             print "\n# self_self_search : Reverse option is set RRRRRRRRRRR ";
                             %fasta_seqs=%{&reverse_sequences( {"$seq_name", "$seq"} )};
                         }else{ %fasta_seqs=("$seq_name", "$seq"); }

                         &do_sequence_search(\%fasta_seqs, "DB=$input_db_file" , "File=$seq_file_name", $create_sso,
                             $single_msp, $over_write, "u=$upper_expect_limit",  "$make_gz_in_sub_dir_opt",
                             "l=$lower_expect_limit", "k=$k_tuple", $make_msp_in_sub_dir_opt, "d=$sub_dir_size");
                         $seq='';
					 }print "\n";
			      }
				  $seq_name=$ori_seq_name;
			  }elsif(eof){
			      $seq.=$_;
				  if($seq=~/\S/ and $seq_name=~/\S/){
					 $seq_file_name="$seq_name\.fa";
					 $seq_file_msp_name="$seq_name\.msp";
                     $first_char=substr("\U$seq_name", 0, $sub_dir_size);
					 if( -s "$first_char\/$seq_file_msp_name" and $over_write !~/o/){
						 print "\n# $first_char\/$seq_file_msp_name already exists ";
					 }else{
                         if($reverse_sequence){ ## reverse the query seqs.
                             %fasta_seqs=%{&reverse_sequences( {"$seq_name", "$seq"} )};
                         }else{ %fasta_seqs=("$seq_name", "$seq"); }
                         &do_sequence_search(\%fasta_seqs, "DB=$input_db_file" , "File=$seq_file_name", $create_sso,
                             $single_msp, $over_write, "u=$upper_expect_limit",  "$make_gz_in_sub_dir_opt",
                             "l=$lower_expect_limit", "k=$k_tuple", $make_msp_in_sub_dir_opt, "d=$sub_dir_size");
                         $seq='';
                     }print "\n";
			      }
			  }elsif(/^(\w+)$/){
				  $seq.=$1;
			  }
		   }
		   close FASTA;
	   }
   }
}

#__________________________________________________________________
# Title     : do_sequence_search
# Usage     : &do_sequence_search("Query_seqs=\%pdb_seq", "DB=$sequence_db_fasta",
#  		         "File=$ARGV[0]", $single_msp, $over_write,
# 	        	 "u=$upper_expect_limit", "l=$lower_expect_limit",
#       		 "k=$k_tuple", $No_processing );
# Function  :
# Example   : &do_sequence_search(\%pdb_seq, $owl_db_fasta, $ARGV[0], $single_msp, $over_write,
#                    "u=$upper_expect_limit", "l=$lower_expect_limit", "k=$k_tuple" );
#
# Keywords  : sequence_search
# Options   :
#             Query_seqs=  for enquiry sequences eg)  "Query_seqs=$ref_of_hash"
#             DB=   for target DB  "DB=$DB_used"
#             File= to get file base(root) name.  "File=$file[0]"
#             m  for MSP format directly from FASTA or Ssearch result than through sso_to_msp to save mem
#             s  for the big single output (msp file output I mean)
#             o  for overwrite existing xxxx.fa files for search
#             c  for create SSO file (sequence search out file)
#             d  for very simple run and saving the result in xxxx.gz format in sub dir starting with one char
#             r  for reverse the query sequence
#             R  for attaching ranges of sequences
#             k= for k-tuple value. default is 1 (ori. FASTA prog. default is 2)
#             u= for $upper_expect_limit
#             l= for $lower_expect_limit
#             a= for choosing either fasta or ssearch algorithm
#             d= for defining the size of subdir made. 2 means it creates
#                    eg, DE while 1 makes D
#             d  for $make_gz_in_sub_dir_opt, putting resultant sso files in gz format and in single char subdir
#             D  for $make_msp_in_sub_dir_opt, convert sso to msp and put in sub dir like /D/, /S/
#             n  for new format to create new msp file format with sso_to_msp routine
#             PVM=  for PVM run of FASTA (FASTA only)
#             M  for machine readable format -m 10 option
#             M= for machine readable format -m 10 option
#             N  for 'NO' do not do any processing but, do the searches only.
#
# Returns   : the names of files created (xxxxx.msp, yyy.msp,,)
# Argument  :
# Version   : 3.9
#----------------------------------------------------------------------------------------
sub do_sequence_search{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	my (@final_out, $add_range, $single_big_msp, $base_name, $create_sso, @nondup,
	   $Single_msp_out_file, %duplicate, $Evalue_thresh, $Score_thresh, @SSO, $sequence_DB,
	   @sso, @temp, $algorithm, $margin, $out_msp_file, @MSP, @final_msp_file_names_out,
	   $upper_expect_limit, $lower_expect_limit, $k_tuple, %seq_input, %MSP, $No_processing,
       $new_format, $PVM_FASTA_run, $over_write, $sub_dir_size, $reverse_sequence );
	my ($E_val) = 5;  ## default 5 <<<<<<<<<<<<<<<<<<<<<

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# DEFAULTS
	#________________________________________
	$k_tuple           =1;  # 1 or 2, 1 is more sensitive
	$algorithm         ='fasta';
	$upper_expect_limit=1;
	$lower_expect_limit=0;
	$Score_thresh      =75;
	$margin            =0;
	$add_range         ='';
	$pwd               =`pwd`; chomp($pwd);
	if($ENV{'PDB40D_FASTA'}){ $sequence_DB  =$ENV{'PDB40D_FASTA'} if $ENV{'PDB40D_FASTA'};
	}else{	print "\n# INFO: Your ENV setting for PDB40D_FASTA doesn't seem to be correct\n";   }

	if($vars{'a'}=~/\S+/){ $algorithm          = $vars{'a'}            };
	if($vars{'u'}=~/\d+/){ $upper_expect_limit = $vars{'u'}            };
	if($vars{'l'}=~/\d+/){ $lower_expect_limit = $vars{'l'}            };
	if($vars{'k'}=~/\d+/){ $k_tuple            = $vars{'k'}            };
	if($vars{'t'}=~/\d+/){ $Score_thresh       = $vars{'t'}            };
	if($vars{'m'}=~/\d+/){ $margin             = $vars{'m'}            };
    if($vars{'d'}=~/\d+/){ $sub_dir_size       = $vars{'d'}            };
	if($vars{'s'}=~/\S+/){ $single_big_msp     = 's'                   };
	if($vars{'DB'}=~/\S+/){ $sequence_DB       = $vars{'DB'}           };
	if($vars{'FILE'}=~/\S+/){ $input_file_name = $vars{'FILE'}; push(@file,$input_file_name) };
	if($vars{'File'}=~/\S+/){ $input_file_name = $vars{'File'}; push(@file,$input_file_name) };
	if($vars{'Query_seqs'}=~/\S+/){ %seq_input = %{$vars{'Query_seqs'}}};
	if($vars{'Query'}=~/\S+/){      %seq_input = %{$vars{'Query'}}};
	if($vars{'u'}    =~/\S+/){ $E_val          = $vars{'u'}            };
	if($vars{'PVM'}  =~/\S+/){ $PVM_FASTA_run  = $vars{'PVM'}; print "\n# PVM opt is set\n";     };
	if($vars{'M'}  =~/\S+/){ $machine_readable = $vars{'M'};           };

	if($char_opt=~/r/){    $reverse_sequence   = 'r' }
	if($char_opt=~/o/){    $over_write         = 'o' }
    if($char_opt=~/c/){    $create_sso         = 'c'; print "# do_sequence_search: \$create_sso is set\n"; }
	if($char_opt=~/s/){    $single_big_msp     = 's'; print "\n# Single file opt is set\n"; }
	if($char_opt=~/m/){    $msp_directly_opt   = 'm' }
	if($char_opt=~/M/){    $machine_readable   = 'M' }
	if($char_opt=~/d/){    $save_in_gz_in_sub_dir  = 'd' }
	if($char_opt=~/D/){$make_msp_in_sub_dir_opt= 'D' } # for simple search and storing msp file
	if($char_opt=~/N/){    $No_processing      = 'N'; $create_sso='c'; }
    if($char_opt=~/R/){    $add_range          = 'R'                   };

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
	# When no %seq is given, but files
	#___________________________________________
	if(@hash==0 and @file > 0){
		print "\n# do_sequence_search: You did not put sequences as in \%seq, but raw sequence file!\n";
		print "        I will run \'open_fasta_files\' sub to fetch sequences to store in \%seq_input\n";
		%seq_input=%{&open_fasta_files(\@file)};
	}else{
		print "\n# do_sequence_search: I will use given seqs in \%seq_input from \%\{\$hash\[0\]\}\n";
		%seq_input=%{$hash[0]};
	}
	my (@list)=keys %seq_input;
	$base_name = ${&get_base_names($input_file_name)};

	print "\n# line:",__LINE__, ", do_sequence_search, \$algorithm => $algorithm, \$base_name:$base_name
	           $input_file_name <--> $sequence_DB\n";

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Controlling which kind of search it should do. Do save_in_gz_in_sub_dir first if d is set
   #_________________________________________________________
   if( $char_opt =~/[dD]/){
	   print "\n# do_sequence_search: You set \'d\' or \'D\' opt with $sequence_DB.\n";
	   my ($seq_name, $seq)= %seq_input;
       my $first_char= substr("\U$seq_name", 0, $sub_dir_size);
	   mkdir ("$first_char", 0777) unless -d $first_char;
	   chdir("$first_char");
	   my $temp_file_name="$seq_name.fa";
	   &write_fasta_seq_by_seq(\%seq_input, $temp_file_name ); ## e makes skip writing when file already
	   $out_file_sso_name="$seq_name\.sso";
	   $out_file_msp_name="$seq_name\.msp";
	   $out_file_gz_name="$seq_name\.msp\.gz";

	   if($char_opt =~/D/){ #### To make MSP file
		  if($machine_readable=~/M/){
             @temp=`$algorithm -m 10 -H  -E $E_val $temp_file_name $sequence_DB $k_tuple`;
		  }else{
			 @temp=`$algorithm -H -E $E_val $temp_file_name $sequence_DB $k_tuple`;
		  }
          print "\n# \@temp has ",scalar(@temp), " lines !\n" if $verbose;

          @msp_hashes_from_temp = @{&open_sso_files(\@temp, $add_range, $create_sso,
		                                            "u=$upper_expect_limit",
		                                            "l=$lower_expect_limit")};
          if(@msp_hashes_from_temp < 1){
              print "\n# do_sequence_search : Error, something is wrong with open_sso_files, LINE=", __LINE__, "\n";
              exit;
          }else{
              print "\n# \@msp_from_temp has ",scalar(@msp_from_temp), " lines !\n";
          }
          @msp_from_temp= values %{$msp_hashes_from_temp[0]};
          print "\n# @msp_from_temp\n" if $verbose;

		  if( !(-s $out_file_gz_name) or $over_write=~/o/){
			  open(MSP, ">$out_file_msp_name");
			  for(@msp_from_temp){    print MSP $_;  }
			  close MSP;
			  if(-s $out_file_gz_name){
			      unlink ($out_file_gz_name);
			      system("gzip $out_file_msp_name"); ## gzipping it
			  }
		  }else{
			  print "\n# Line No. ", __LINE__,", $out_file_gz_name already exists  (do_sequence_search)\n";
		  }
	   }else{ ### To make gzipped SSO files
		  system(" $algorithm -m 10 -H  -E $E_val $temp_file_name $sequence_DB $k_tuple > $out_file_sso_name");
		  system("gzip $out_file_sso_name");
	   }
	   unlink("$seq_name.fa");
	   print "\n# Sub dir $first_char has been made, finishing do_sequence_search\n";
	   chdir ('..');
	   goto EXIT;
   }


   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # This is the big single MSP output
   #______________________________________________
   $Single_msp_out_file="$base_name\.msp" if($single_big_msp eq 's');
   if(-s $Single_msp_out_file and $char_opt !~/o/){
	   print "\n# $Single_msp_out_file exists, skipping \n";
	   push(@final_out, $Single_msp_out_file);
	   return(\@final_out);
   }else{
	   $char_opt .='o';
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Check if it is necessary to write sequences.fa files
   #______________________________________________________
   if($char_opt=~/o/){
	  &write_fasta_seq_by_seq(\%seq_input); ## e makes skip writing when file already
   }else{
	  &write_fasta_seq_by_seq(\%seq_input, 'e'); ## e makes skip writing when file already
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
   #  When, you didn't use "DB=$XXX" and "File=$FXXX" format, first file input is DB etc
   #_______________________________________________________________________________________

   if($sequence_DB=~/^$/){
	  print "\n# FATAL: do_sequence_search: You did not use \"DB=\$XXX\" format\n"; exit   };

   print "\n# Finished writing the enquiry fasta files from \%seq_input by write_fasta_seq_by_seq";
   print "\n# I am in do_sequence_search sub, Target database used :  $sequence_DB with seqs of \'@list\'\n";


   for($j=0; $j< @list; $j++){  # @list has sequence names
	   my @temp;
	   my $each_seq_fasta="$list[$j]\.fa";
	   unless(-s $each_seq_fasta){   print "\n# do_sequence_search: $each_seq_fasta does not exist, error\n"; exit }
	   print "\n# Found $each_seq_fasta is searched against $sequence_DB\n";
	   $out_msp_file="$list[$j]\.msp";
	   $out_sso_file ="$list[$j]\.sso";
	   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	   # If files already exist
	   #__________________________________________
	   if( -s $out_msp_file and $char_opt !~/o/ ){
		   print "\n# File: $out_msp_file exists, skipping, to overwrite use \'o\' opt";

		   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		   # c opt for creating SSO file
		   #____________________________________
		   if($create_sso=~/c/){
		      if($machine_readable=~/M/){
				 @temp=`$algorithm -m 10 -H  -E $E_val $each_seq_fasta $sequence_DB $k_tuple`;
			  }else{
				 @temp=`$algorithm -H -E $E_val $each_seq_fasta $sequence_DB $k_tuple`;
			  }
			  if(@temp < 20){
				  print "\n# OUTPUT of fasta is too small, error \n"; print chr(7);
				  exit;
			  }
			  open(SSO, ">$out_sso_file");
			  print SSO @temp;
			  print "\n# $out_sso_file is created";
			  close SSO;
		   }
		   push(@final_out, $out_msp_file);
		   unlink($each_seq_fasta);
	   }
	   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	   # If files DONT  exist
	   #__________________________________________
	   else{  ## -E is for e value cutoff. -b is for num of seq fetched
		   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		   #  K-tuple is  set  to  1 by default
		   #
		   #  If xxxx.sso exists, let's skip running fasta or ssearch
		   #____________________________________________________________

		   if(-s $out_sso_file and $char_opt !~/o/ ){
			   open(SSO_ALREADY, "$out_sso_file");
			   @temp=<SSO_ALREADY>;
			   close(SSO_ALREADY);
		   }else{
			   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			   #  NORMAL Default FASTA run comes here
			   #________________________________________
			   if($machine_readable=~/M/){   print "\n# do_sequence_search: You put \'M\' opt \n";
			      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			      # PVM FASTA run
			      #__________________________________________
				  if($PVM_FASTA_run=~/PVM/){
				     unless(-s $ENV{'PVM_HOSTFILE'}){
				         print "\n# do_sequence_search: $ENV{'PVM_HOSTFILE'} does not exsists error\n";
						 exit;
				     }
					 open(PVM_CSH, ">pvm_fasta.csh"); print "\n# pvm_fasta.csh is created \n";
				     print PVM_CSH "\#\!\/bin\/csh\n";
					 print PVM_CSH "pvm $ENV{'PVM_HOSTFILE'} \<\< \'eof\'\n\'eof\'\n";
					 print PVM_CSH '/gn0/jong/App/Fasta/pvcompfa -m 10 -H -E ', " $E_val ",
									" $each_seq_fasta ", "$sequence_DB $k_tuple \> temp.sso\n";

					 print PVM_CSH "\npvm \<\< \'eof\'\n";
					 print PVM_CSH "halt\n\'eof\'\n";
					 close PVM_CSH;
					 system(" csh pvm_fasta.csh");
					 open(TEMP_SSO, "temp.sso");
					 @temp=<TEMP_SSO>;
					 close TEMP_SSO;
				  }else{
					 @temp=`$algorithm -m 10 -H  -E $E_val $each_seq_fasta $sequence_DB $k_tuple`;
				  }
			   }else{
			      if($PVM_FASTA_run=~/PVM/){
				     unless(-s $ENV{'PVM_HOSTFILE'}){
				         print "\n# do_sequence_search: $ENV{'PVM_HOSTFILE'} does not exsists error\n";
						 exit;
				     }
					 open(PVM_CSH, ">pvm_fasta.csh");
				     print PVM_CSH "\#\!\/bin\/csh\n";
					 print PVM_CSH "pvm $ENV{'PVM_HOSTFILE'} \<\< \'eof\'\n\'eof\'\n";
					 print PVM_CSH '/gn0/jong/App/Fasta/pvcompfa -H -E ', " $E_val ",
									" $each_seq_fasta ", "$sequence_DB $k_tuple \> temp.sso \n";
					 print PVM_CSH "\npvm \<\< \'eof\'\n";
					 print PVM_CSH "halt\n\'eof\'\n";
					 close PVM_CSH;
					 system("csh pvm_fasta.csh");
					 open(TEMP_SSO, "temp.sso");
					 @temp=<TEMP_SSO>;
					 close TEMP_SSO;
			      }else{
			         @temp=`$algorithm -H -E $E_val $each_seq_fasta $sequence_DB $k_tuple`;
				  }
			   }

			   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			   # c opt for creating SSO file
			   #____________________________________
			   if($create_sso=~/c/){
				  if(@temp < 20){
					  print "\n# OUTPUT of fasta is too small, error \n"; print chr(7);
					  exit;
				  }
				  open(SSO, ">$out_sso_file");
				  print SSO @temp;
				  print "\n# $out_sso_file   is created because of \"c\" or \"N\" option you set ";
				  close SSO;
			   }
		   }

		   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
		   # If $char_opt is m
		   #____________________________________________________________
		   if($char_opt=~/m/){  # make msp files directly from search not going through sso_to_msp
			   my @msp_hashes_from_temp = @{&open_sso_files(\@temp, $add_range,
			                                                "u=$upper_expect_limit",
			                                                "l=$lower_expect_limit")};
			   my @msp_from_temp= values %{$msp_hashes_from_temp[0]};
			   $MSP{$out_msp_file} = \@msp_from_temp;
			   unlink($each_seq_fasta);
			   next;
		   }elsif($No_processing !~/N/){ ## When sso output is not directly converted to MSP
				 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				 # REmoving some junk lines  in SSO output
				 #_________________________________
				 for($o=0; $o < @temp; $o++){
					 if($temp[$o]=~/^[ \w\-]+/){                    splice(@temp, $o, 1); $o--;
					 }elsif($temp[$o]=~/^\; mp/){                   splice(@temp, $o, 1); $o--;
					 }elsif($temp[$o]=~/^\; pg/){                   splice(@temp, $o, 1); $o--;
					 }elsif($temp[$o]=~/^\; fa_[ozi]/){             splice(@temp, $o, 1); $o--;
					 }elsif($temp[$o]=~/^\; sq_type$/){             splice(@temp, $o, 1); $o--;
					 }
				 }
				 unlink($each_seq_fasta);
				 if(@temp < 20){
					 print "\n# FATAL: OUTPUT of FASTA is too small (less than 20 byte), error\n";
					 print "\n# @temp\n";
					 print chr(7);
					 exit;
				 }
				 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				 # opt s is for the big single output
				 #___________________________________
				 if($single_big_msp eq 's'){  push(@SSO, \@temp); ## <<------------ to go to (3)
				 }else{
					 $msp_out_file="$list[$j]\.msp";
					 push(@final_out, @{&sso_to_msp(\@temp, "l=$lower_expect_limit",
										 "u=$upper_expect_limit", $msp_out_file,
										 $create_sso, $add_range, "t=$Score_thresh",
										 "e=$Evalue_thresh", "m=$margin", $new_format )} );
					 if($char_opt=~/c/){ ##  create SSO file (sequence search out file
						$out_sso_file ="$list[$j]\.sso";
						open(SSO, ">$out_sso_file");
						print SSO @temp;
						print "\n# $out_sso_file is created \n";
						close SSO;
					 }
				 }
		   }else{   # endof if($char_opt=~/m/){ }
				 print "\n# do_sequence_search: You set \'N\' option for NO processing of the results\n";
		   }
	   }
   } # end of for($j=0; $j< @list; $j++){

   if($char_opt=~/m/){  # make msp files directly from search not going through sso_to_msp
	   if($single_big_msp=~/s/){
		  open(SINGLE_BIG_MSP, ">$Single_msp_out_file");
		  @MSP= keys %MSP;
		  for($m=0; $m< @MSP; $m++){
			 print SINGLE_BIG_MSP @{$MSP{$MSP[$m]}}, "\n";
		  }
		  close(SINGLE_BIG_MSP);
		  push(@final_msp_file_names_out, $Single_msp_out_file);
		  return(\@final_msp_file_names_out);
	   }else{
		  @MSP= keys %MSP;
	      for($t=0; $t <  @MSP; $t++){
			 open(SING_MSP, ">$MSP[$t]");
	         print SING_MSP @{$MSP{$MSP[$t]}}, "\n";
			 close(SING_MSP);
			 push(@final_msp_file_names_out, $Single_msp_out_file);
		  }
		  return(\@final_msp_file_names_out);
	   }
   }else{
	   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	   # (3) This now processes the @SSO which has all the sso Single SSO
	   #________________________________________________________
	   if($single_big_msp =~ /s/){
		  push(@final_out, @{&sso_to_msp(@SSO, $Single_msp_out_file, $create_sso,
			 "u=$upper_expect_limit", "l=$lower_expect_limit",
			 $add_range, "m=$margin", $new_format)} );
		  if($char_opt=~/c/){ ##  create SSO file (sequence search out file
			 $out_sso_file ="$base_name\.sso";
			 open(SSO, ">$out_sso_file");
			 for($i=0; $i< @SSO; $i++){
				print SSO @{$SSO[$i]}, "\n";
			 }
			 print "\n# $out_sso_file is created \n";
		  }
		  close(SSO);
	   }
	   @nondup = grep { ! $duplicate{$_}++ } @final_out;
	   return(\@nondup);
   }
   EXIT:
}








#________________________________________________________________________
# Title     : assign_options_to_variables
# Usage     : &assign_options_to_variables(\$input_line);
# Function  : Assigns the values set in head box to the variables used in
#             the programs according to the values given at prompt.
#             This produces global values.
#             When numbers are given at prompt, they go to @num_opt
#              global variable. %vars global option will be made
#
# Example   : When you want to set 'a' char to a variable called '$dummy' in
#             the program, you put a head box commented line
#             '#  $dummy    becomes  a  by  -a '
#             Then, the parse_arguments and this sub routine will read the head
#             box and assigns 'a' to $dummy IF you put an argument of '-a' in
#             the prompt.
# Warning   : This is a global vars generator!!!
# Keywords  :
# Options   : '#' at prompt will make a var  $debug set to 1. This is to
#              print out all the print lines to make debugging easier.
# Returns   : Some globaly used variables according to prompt options.
#             @num_opt,
#
# Argument  : None.
# Version   : 2.6
#--------------------------------------------------------------------
sub assign_options_to_variables{
  my($i, %vars, $j, $op, $z, $n, $symb, $value, $var, %val, @val, $ARG_REG,
	 $option_table_example, @input_options, $first_border_and_title, $sym, @arg);

  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #      Defining small variables for option table reading
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  my($g)='gets';                my($if)='if';
  my($is)='is';                 my(@input_files);
  my($o)='or';   my(@arguments) = sort @ARGV;

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  Assigning global arguments(@num_opt, %vars) variables
  #_______________________________________________________________
  for($i=0; $i< @arguments; $i++){
	 if(($arguments[$i]=~/^(\-?\d+[\.\d+]?)$/)&&   ### it mustn't be a file
		( !(-f $arguments[$i]) ) ){                ### getting NUM opt
		push(@num_opt, $1);
	 }elsif( $arguments[$i]=~/^(\S+)=(\S+)$/){
		$vars{$1}=$2;
	 }
  }

  #""""""""""""""""""""""""""""""""""""""""""""""""""
  #   Some DEFAULT $debug variables for debugging purposes
  #""""""""""""""""""""""""""""""""""""""""""""""""""
  &set_debug_option;

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #   The main processing of self
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  open(SELF, "$0");    ## opens the program you ran to get the options table.
  while(<SELF>){

	  if( $first_border_and_title > 6 ){  ## This is to make it read only the first headbox.
		  last;                            #  $first_border_and_title is an incremental counter.
	  }elsif( /^ *#[_\*\-]{15,}$/ and /^ *# *[Tt][itle]*[ :]*/ ){
		  $first_border_and_title++;
		  print __LINE__, "# assign_options_to_variables : Title line found\n" if $debug eq 1;
	  }elsif(/^ {0,5}# {1,50}[\$\%\@].+$/){
		  $op = $&;  ## $op is for the whole input option line which has $xxxx, @xxx, %xxxx format
		  $op =~ s/^( *\# *)(\W\w+.+)$/$2/;  ## This is removing '#  ' in the line.
		  $op =~ s/^(\W\w+.+)(\s+\#.*)$/$1/;  ## This is removing any comments in the line.
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 ## matching the following line input format.
			 ## $av_sc_segment     becomes    a  by  a  # To smooth the SC rates. Gets the averages of
			 ## $ARG_REG is for arguments regular expression variable.
			 ##  This reg. exp. matches = 'a or A or E or e' part
			 ##  which represents alternative prompt arguments possibilities. \=$b$g$is$e$set
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 $ARG_REG ='(\S*) *[or=\,]* *(\S*) *[=or\,]* *(\S*) *[=or\,]* *(\S*) *[=or\,]* *(\S*)';
			 if($op=~/^([\$\@\%])([\w\-]+) {0,20}[=|$g|$is] *[\$\@\%]*([\- \w\.\d]+) *[bB]y +$ARG_REG/){
							 ## $sym     $var        becomes          a [$a...]       by       a -a -A
				  my $sym = $1;  #### The symbols like ($, @, %), '$' in the above.
				  my $var = $2;  #### Actual variable name 'var' from $var, 'av_sc_segment' in the above.
				  my $val = $3;  #### The becoming value  first 'a' in the above.
				  my @arg = ($4, $5, $6, $7, $8);  ## The alternative prompt arguments, second 'a' in the above..
			      print "\n $sym $var $val \n" if $debug==1;
			      print "\n \@arg are @arg \n" if $debug==1;

				  #""""""""""""""""""""""""""""""""""""""""""""""""""""
				  #  Going through the PROMPT args.
				  #""""""""""""""""""""""""""""""""""""""""""""""""""""
				  for($z=0; $z < @arguments; $z++){     ## $arguments[$z]  is from @ARGV
					  if($arguments[$z]=~/^\-\w+$/){
						  $arguments[$z] =~ s/\-//;
					  }
					  for ($i=0; $i < @arg; $i ++ ){
						 if( ("$arg[$i]" eq "$arguments[$z]" )&& ($arg[$i] !~ /\=/)
							 && ($sym eq '$') ){
							 ${"$var"}="$val";
							 if($debug == 1){
								 print __LINE__," \$${var} is set to \"$1\"\n";
							 }

						 }#'''''''''''''''' $arg = by s=  syntax ~~~~~~~~~~~~~~~~~~~~~~~~~~~
						 elsif( ( $arg[$i] =~ /^(\w+) *\=/ ) &&
							( $arguments[$z] =~ /^${1}= *([\w\.*\-*]+)$/) &&
							( $sym eq '$') ){
							  ${"$var"}="$1";
							  if($debug eq 1){ print __LINE__,"\$${var} is set to \"$1\"\n";  }
						 }
					  }
				  }
			  }
		}
	}
}
#________________________________________________________________________
# Title     : default_help
# Usage     : &default_help2;  usually with 'parse_arguments' sub.
# Function  : Prints usage information and others when invoked. You need to have
#             sections like this explanation box in your perl code. When invoked,
#             default_help routine reads the running perl code (SELF READING) and
#             displays what you have typed in this box.
#             After one entry names like # Function :, the following lines without
#             entry name (like this very line) are attached to the previous entry.
#             In this example, to # Function : entry.
# Example   : &default_help2; &default_help2(\$arg_num_limit);   &default_help2( '3' );
#             1 scalar digit for the minimum number of arg (optional),
#             or its ref. If this defined, it will produce exit the program
#             telling the minimum arguments.
# Warning   : this uses format and references
# Keywords  :
# Options   :
# Returns   : formated information
# Argument  :
# Version   : 3.3
#--------------------------------------------------------------------
sub default_help{
  my($i, $perl_dir, $arg_num_limit, $head ,$arg_num_limit,
	  @entries, @entries_I_want_write );
  my($logname)=getlogin();
  my($pwd)=`pwd`;
  my($date)=`date`;
  chomp($date,$pwd);
  my($not_provided)="--- not provided ---\n";
  my($file_to_read) = $0;

  for($i=0; $i < @_; $i ++){
	  if((ref($_[$i]) eq 'SCALAR')&&(${$_[$i]} =~ /^\d$/)){
		  $arg_num_limit = ${$_[$i]};  }
	  elsif( (!(ref($_[$i]))) && ($_[$i] =~ /^\d$/)){
		  $arg_num_limit = $_[$i];     }
  }
  my %entries = %{&read_head_box(\$file_to_read )};
  if($option_tb_found ==1){
	 @option_tb=@{&read_option_table(\$file_to_read)};
  }

  @entries = keys %entries;
  foreach $help_item (@entries){
	  ${$help_item}= $not_provided if( ${$help_item}=~/^[\W]*$/  and  !defined(${$help_item}) );
  }
  #""""""""""""""""""""""""""""""""""""""""
  #########  Writing the format <<<<<<<<<<<
  #""""""""""""""""""""""""""""""""""""""""
  $~ =HEADER_HELP;
  write;   ## <<--  $~ is the selection operator
  $~ =DEFAULT_HELP_FORM;

  @entries_I_want_write=sort keys %entries;

  for( @entries_I_want_write ){  write  }

  print chr(7);  print "_"x72,"\n\n";

  if(@ARGV < $arg_num_limit){ print "\* $0 fataly needs $arg_num_limit arguments\n\n" }

  if(  $option_tb_found == 1){
	 #########  Printing the OPTION table contents <<<<<<<<<<<<
	 print "  Press \"Return\" key to see what options $logname\'s \n\n    \"$0\" take... \n";
		 $key_press=getc();
	 print @option_tb, "\n"x2 if(@option_tb > 0);
  }
format HEADER_HELP  =
_____________________________________________________________________
		  __  __      ______     __          _____
		 /\ \/\ \    /\  ___\   /\ \        /\  _ `\
		 \ \ \_\ \   \ \ \__/   \ \ \       \ \ \L\ \
		  \ \  _  \   \ \  _\    \ \ \       \ \ ,__/
		   \ \ \ \ \   \ \ \/___  \ \ \_____  \ \ \/
		    \ \_\ \_\   \ \_____\  \ \______\  \ \_\
		     \/_/\/_/    \/_____/   \/______/   \/_/ V 3.1`
_____________________________________________________________________
.
format DEFAULT_HELP_FORM =
 @<<<<<<<<<: @*
 $_        $entries{$_}
.
}
#________________________________________________________________________
# Title     : set_debug_option
# Usage     : &set_debug_option;
# Function  : If you put '#' or  '##' at the prompt of any program which uses
#             this sub you will get verbose printouts for the program if the program
#             has a lot of comments.
# Example   : set_debug_option #    <-- at prompt.
# Warning   :
# Keywords  :
# Options   : #   for 1st level of verbose printouts
#             ##  for even more verbose printouts
# $debug  becomes 1 by '#'  or '_'
# $debug2 becomes 1 by '##'  or '__'
#
# Returns   :  $debug
# Argument  :
# Version   : 1.8
#--------------------------------------------------------------------
sub set_debug_option{
  my($j, $i, $level);
  unless( defined($debug) ){
	 for($j=0; $j < @ARGV; $j ++){
		 if( $ARGV[$j] =~/^(_+)$|^(#+)$/){ # in bash, '#' is a special var, so use '_'
			 print __LINE__," >>>>>>> Debug option is set by $1 <<<<<<<<<\n";
			 $debug=1;
				  print chr(7);
			 print __LINE__," \$debug  is set to ", $debug, "\n";
			 splice(@ARGV,$j,1); $j-- ;
			 $level = length($1)+1;
			 for($i=0; $i < $level; $i++){
				 ${"debug$i"}=1;
				 print __LINE__," \$debug${i} is set to ", ${"debug$i"}, "\n";
			 }
		 }
	 }
  }
}
#___________________________________________________________
# Title     : get_seq_fragments
# Usage     : @seq_frag=&get_seq_fragments(\%msf, @RANGE);
# Function  : gets sequence(string) segments with defined
#             ranges.
# Example   :
#  %test=('seq1', '1234AAAAAAAAAAAaaaaa', 'seq2', '1234BBBBBBB');
#  @range = ('1-4', '5-8');
#
#  %out = %{&get_seq_fragments(\%test, \@range)};
#  %out => (seq1_5-8   AAAAA
#           seq2_5-8   BBBBB
#           seq1_1-4    1234
#           seq2_1-4    1234 )
#
# Warning   :
# Keywords  : get_sequence_fragments,
# Options   : _  for debugging.
#             #  for debugging.
#             l=  for min seqlet length
#             r  for adding ranges in the seq names
#
# Returns   :
# Argument  :
# Version   : 1.8
#-------------------------------------------------------
sub get_seq_fragments{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
   my $min_seqlet_size=10;
   my @vars=keys %vars;
   my $no_range_in_name=1;
   for($i=0; $i< @vars; $i++){
	   if($vars[$i] eq 'l'){
		  $min_seqlet_size=$vars{$vars[$i]};
	   }
   }
   if($char_opt=~/v/){ print "\n \$char_opt is $char_opt  @char_opt\n"; }
   if($char_opt=~/n/){ $no_range_in_name = 1 }
   if($char_opt=~/r/){ $no_range_in_name = 0 }

   print "\nget_seq_fragments \$no_range_in_name is $no_range_in_name \n";
   for($i=0; $i< @hash; $i++){
	 my (%out_frag, $frag_name, $range_start, $range_end, @out_hash);
	 my %seqs = %{$hash[$i]};
	 my @names = keys %seqs;
	 if(@names==1){
	    for($j=0; $j < @names; $j++){
		   my $seq_name = $names[$j];
		   my $seq = $seqs{$seq_name};
		   for($k=0; $k< @range; $k++){
			  my $range = $range[$k];
			  if($no_range_in_name==1){
				 $frag_name = "$seq_name";
			  }else{
			     $frag_name = "$seq_name\_$range";
			  }
			  #if(length($frag_name)>14 ){
			  #	 $frag_name ='x'."${j}_${range}";
		      #}
			  ($range_start, $range_end)=$range=~/(\d+\.?\d*)\-(\d+\.?\d*)/;
			  my $frag_len = $range_end-$range_start+1;
			  if($frag_len < $min_seqlet_size){
			     next;
			  }
			  my $fragment = substr($seq, $range_start-1, $frag_len);
			  $out_frag{$frag_name}=$fragment;
		   }
		}
		push(@out_hash,  \%out_frag);
	 }elsif(@names > 1){
	    for($k=0; $k< @range; $k++){
		  my %out_frag=();
	      my $range=$range[$k];
		  ($range_start, $range_end)=$range=~/(\d+\.?\d*)\-(\d+\.?\d*)/;
	      my $frag_len = $range_end-$range_start+1;
		  if($frag_len < $min_seqlet_size){
		     next;
		  }
	      for($j=0; $j < @names; $j++){
	         my $seq_name=$names[$j];
			 my $seq = $seqs{$seq_name};
		     if($no_range_in_name==1){
				 $frag_name = "$seq_name";
			 }else{
			     $frag_name = "$seq_name\_$range";
			 }
			 #if(length($frag_name)>15 ){
			 #	$frag_name ='x'."${j}_${range}";
		     #}
			 if($range_start==0){ $range_start++; } ## This is a bugfix
			 my $fragment = substr($seq, $range_start-1, $frag_len);
			 $out_frag{$frag_name}=$fragment;
		  }
		  push(@out_hash, \%out_frag);
		}
	 }
   }
   if(@out_hash > 1){ return(@out_hash)
   }elsif(@out_hash==1){ return($out_hash[0]) }
}
#__________________________________________________________________________
# Title     : if_file_older_than_x_days
# Usage     : if( ${&if_file_older_than_x_days($ARGV[0], $days)} > 0){
# Function  : checks the date of last modi of file given and compares with
#             present time. Substracts diff and returns the actual diff days.
# Example   :
# Keywords  : how_old_file, how_old, is_file_older_than_x_days, file_age,
#             file_age_in_days,
# Options   :
# Returns   : the actual days older, so NON-ZERO, otherwise, 0
# Version   : 1.3
#----------------------------------------------------------------------------
sub if_file_older_than_x_days{
	if(@_ < 2){ print "\n# FATAL: if_file_older_than_x_days needs 2 args\n"; exit; }
	my $file=${$_[0]} || $_[0];
	my $days=${$_[1]} || $_[1];
	my ($new_idx_file, $how_old_days);
	unless(-s $file){
	    print "\n# FATAL, nearly!: if_file_older_than_x_days: $file does NOT exist !\n";
		$new_idx_file=${&make_seq_index_file($file)};
		print "        if_file_older_than_x_days called make_seq_index_file to make $new_idx_file\n";
	}

	$how_old_days=(localtime(time- (stat($file))[9]))[3];
	if($how_old_days > $days){
		print "\n# if_file_older_than_x_days: $file is older than $days\n";
		return(\$days);
	}else{
		print "\n# if_file_older_than_x_days: $file is NOT older than $days\n";
		return(0);
	}
}
#________________________________________________________________________
# Title     : parse_arguments
# Usage     : &parse_arguments; or  (file1, file2)=@{&parse_arguments};
# Function  : Parse and assign any types of arguments on prompt in UNIX to
#             the various variables inside of the running program.
#             This is more visual than getopt and easier.
#             just change the option table_example below for your own variable
#             setttings. This program reads itself and parse the arguments
#             according to the setting you made in this subroutine or
#             option table in anywhere in the program.
# Example   : &parse_arguments(1);
#             @files=@{&parse_arguments(1)};
# Warning   : HASH and ARRAY mustn't be like = (1, 2,3) or (1,2 ,3)
# Keywords  :
# Options   : '0'  to specify that there is no argument to sub, use
#              &parse_arguments(0);
#             parse_arguments itself does not have any specific option.
#             '#' at prompt will make a var  $debug set to 1. This is to
#              print out all the print lines to make debugging easier.
#             'e=xxxx' for filtering input files by extension xxxx
#
# Returns   : Filenames in a reference of array
#             and input files in an array (file1, file2)=@{&parse_arguments};
# Argument  : uses @ARGV
# Version   : 1.8
#--------------------------------------------------------------------
sub parse_arguments{
  my( $c, $d, $f, $arg_num, $option_table_seen, $n, $option_filtered,
		$option_table_example, $input_line, @input_files,
		$extension);
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #   Checks if there were arguments
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if( @ARGV < 1 ){ #<-- If Argument is not given at prompt
	  for(@_){
		 if($_ eq '0'){
			 last;
		 }else{
			 print "\n \"$0\" requires at least one Argument, suiciding.\n\n";
			 print chr(7); #<-- This is beeping
			 print "  To get help type \"$0  h\"\n\n\n ";
			 exit;
		 }
	  }
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  Checking some input options like 'e=txt' for extension filtering
  #_____________________________________________________________________
  for($i=0; $i< @_; $i++){
	  if($_[$i]=~/e=(\S+)/){
		  push(@extension, $1);
	  }
  }

  #""""""""""""""""""""""""""""""""""""""""""""""""""
  #   Some DEFAULT $debug variables for debugging purposes
  #""""""""""""""""""""""""""""""""""""""""""""""""""
  &set_debug_option;
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  If there is only one prompt arg. and is [-]*[hH][elp]*, it calls
  #   &default_help and exits
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if( ( @ARGV == 1 ) && ($ARGV[0] =~ /^[\-]*[hH\?][elp ]*$/) ){
		&default_help;
		exit;
  }
  for($f=0; $f < @ARGV; $f++){
	 if( $ARGV[$f] =~ /\w+[\-\.\w]+$/ and -f $ARGV[$f] ){
		 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		 # When extension is defined, filter files by it
		 #____________________________________________
		 if(@extension > 0){
		     for($e=0; $e < @extension; $e++){
				 $extension=$extension[$e];
				 if($ARGV[$f]=~/\S\.$extension/){
					 push(@input_files, $ARGV[$f] );
				 }else{ next }
			 }
		 }else{
			 push(@input_files, $ARGV[$f] );
			 next;
		 }
	 }
  }

  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #     Reading the running program script
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  &assign_options_to_variables;
  if($HELP == 1){ &default_help }
  return(\@input_files);
}
#________________________________________________________________________
# Title     : read_option_table
# Usage     :
# Function  : Reads the option table made by Jong in any perl script. The
#             option table is a box with separators.
# Example   :
# Warning   :
# Keywords  :
# Options   :
# Returns   :
# Argument  :
# Version   : 1.0
#--------------------------------------------------------------------
sub read_option_table{
	my($table_found, @option_tb, $head);
	 open(SELF, "${$_[0]}");
	 while(<SELF>){
		if( (/^ *#+/) && ( $table_found== 1) ){
		  push (@option_tb, "$_");
		}elsif( ($table_found != 1)&&(/^ *\#+ *[Oo]ption *[Tt]able */) ){
			$table_found=1; $head="############## Option Table for $logname\'s \"$0\"\n"; ##
			push(@option_tb, $head);
		}
		if( ($table_found==1)&&(/^ *###################+ *$/)){  ### to find the end point of reading
			$table_found =0; last; }
	 }
	 return(\@option_tb);
}
#______________________________________________________________________________________
# Title     : open_sso_files
# Usage     :  @sso=@{&open_sso_files(@file, $add_range, $add_range2, "u=$upper_expect_limit",
#			                            "l=$lower_expect_limit", "m=$margin", $new_format)};
# Function  : This reads the parseable( -m 10 option)
#              and non-parseable form of ssearch program output
#             If you give 5 files, it produces 5 hashes as a ref of array.
#             This understands xxxx.gz files.
#             This reads FASTA -m 10 output, too.
# Example   :
#  717    0         0.343  16    373    EC1260_16-373              74    434    YBL6_YEAST_74-434
#  348    9e-16     0.500  113   233    EC1260_113-233             27    146    YDBG_ECOLI_27-146
#  472    2.9e-08   0.271  13    407    EC1260_13-407              148   567    YHJ9_YEAST_148-567
#  459    1.9e-22   0.260  1     407    EC1260_1-407               65    477    YLQ6_CAEEL_65-477
#  452    4.5e-14   0.275  1     407    EC1260_1-407               103   537    YSCPUT2_103-537
#  1131   0         0.433  1     407    EC1260_1-407               112   519    ZMU43082_112-519
#
# Warning   : By default, the SW score comes to the first
#             If expect value is not found, it becomes '0'
#             By default, the offset of seq match with a seq name like seq_30-40
#               will be 30 not 1.
#             It ignores special chars like , : .prot in the name (eg, AADF_FASDF: will be AADF_FASDF)
# Keywords  : open_ssearch_output_files, ssearch_output, ssearch, FASTA,
# Options   : _  for debugging.
#             #  for debugging.
#             u= for upper E value limit
#             l= for lower E value limit
#             r  for attaching ranges to out seq names (eg> HI0001_1-20 as a key)
#             U  for making the matched seqname to upppercase
#             L  for making the matched seqname to lowercase
#             R  for attaching ranges to out seq names for both TARGET and MATCH
#             n  for new format (msp2)
#             a  for getting alignments of the pair
#
# Version   : 4.2
# Enclosed  :
#
#   >>MG032 ATP-dependent nuclease (addA) {Bacillus subtilis  (666 aa)
#    Z-score: 88.3 expect()  1.9
#   Smith-Waterman score: 77;  27.143% identity in 70 aa overlap
#
#           30        40        50        60        70        80
#   MJ0497 RSAGSKGVDLIAGRKGEVLIFECKTSSKTKFYINKEDIEKLISFSEIFGGKPYLAIKFNG
#                                        : .. ...  . .:.:::. :: : ..:
#   MG032  HDKVRYAFEVKFNIALVLSINKSNVDFDFDFILKTDNFSDIENFNEIFNRKPALQFRFYT
#        200       210       220       230       240       250
#
#           90       100             110       120       130
#   MJ0497 EMLFINPFLLSTNGK------NYVIDERIKAIAIDFYEVIGRGKQLKIDDLI
#          .   ::   :: ::.      : ....... . ::. . :
#   MG032  K---INVHKLSFNGSDSTYIANILLQDQFNLLEIDLNKSIYALDLENAKERFDKEFVQPL
#        260          270       280       290       300       310
#
# Parseable form -m 10 option =========================================
#   >>>MJ0497.fa, 133 aa vs GMG.fa library
#   ; pg_name: Smith-Waterman (PGopt)
#   ; pg_ver: 3.0 June, 1996
#   ; pg_matrix: BL50
#   ; pg_gap-pen: -12 -2
#   >>MG032 ATP-dependent nuclease (addA) {Bacillus subtilis
#   ; sw_score:  77
#   ; sw_z-score: 88.3
#   ; sw_expect    1.9
#   ; sw_ident: 0.271
#   ; sw_overlap: 70
#   >MJ0497 ..
#   ; sq_len: 133
#   ; sq_type: p
#   ; al_start: 58
#   ; al_stop: 121
#   ; al_display_start: 28
#----------------------------------------------------------------------------
sub open_sso_files{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	my (@out_refs, @SSO, $create_sso, $parseable, @OUT, @temp_sso_lines,
		%match, $attach_range_in_names, $margin, $uppercase_seq_name,
		$lowercase_seq_name, $target_seq, $new_format, $get_alignment,
		$pvm_version_fasta_out, $original_target_seq);

	my ($upper_expect_limit, $lower_expect_limit)=(50,0);

	if($char_opt=~/R/){  $attach_range_in_names2=1; };
	if($char_opt=~/r2/){ $attach_range_in_names =1; $attach_range_in_names2=1 };
	if($char_opt=~/r/){  $attach_range_in_names =1; };
	if($char_opt=~/c/){  $create_sso   ='c' ; print "\n# open_sso_files: \$create_sso is set";};
	if($char_opt=~/n/){  $new_format   ='n' };
	if($char_opt=~/a/){  $get_alignment='a' };
	if($char_opt=~/U[pperPPER]*/){ $uppercase_seq_name='U' };
	if($char_opt=~/L[owerOWER]*/){ $lowercase_seq_name='L' };
	if($vars{'u'}=~/([\.\d]+)/){ $upper_expect_limit = $vars{'u'} };
	if($vars{'l'}=~/([\.\d]+)/){ $lower_expect_limit = $vars{'l'} };
	if($vars{'m'}=~/\d+/){ $margin = $vars{'m'} };

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# opening file input (can handle .gz  files)
	#_______________________________________________
    if(@file < 1 and @array > 0){
         for($i=0; $i< @array; $i++){
              @sso=@{$array[$i]};
         }
         print "\n# \@sso has ", scalar(@sso), " lines. \n" if $verbose;
         if(@sso > 3000){ # if @sso is very big, I remove the useless contents
             print "\n# open_sso_files: size of \@sso for $file[$i] exceeds 3000 lines, ", scalar(@sso), " !!! \n";
         }
         push(@OUT, &read_sso_lines(\@sso, $create_sso, $attach_range_in_names, $attach_range_in_names2,
                         $new_format, $get_alignment) );
    }else{
         for($i=0; $i< @file; $i++){
              if($file[$i]=~/\S+\.\gz$/ or -B $file[$i]){  ## if file has xxxx.gz extension
                  my (@sso);
                  @sso=`gunzip -c $file[$i]`;
                  if(@sso < 30){  @sso=`zcat $file[$i]`; }      # if zcat fails to produce output use gunzip -c
                  if(@sso > 3000){ # if @sso is very big, I remove the useless contents
                      print "\n# open_sso_files: size of \@sso for $file[$i] exceeds 3000 lines, ", scalar(@sso), " !!! \n";
                  }
                  push(@OUT, &read_sso_lines(\@sso, $create_sso, $attach_range_in_names, $attach_range_in_names2,
                                  $new_format, $get_alignment) );
              }else{
                  print "\n# openning text file format xxxx.sso $file[$i]";
                  open(SSO, "$file[$i]") or die "\n# open_sso_files: Failed to open $file[$i]\n";
                  my @sso=<SSO>;
                  if(@sso < 30){  @sso=`zcat $file[$i]`; }      # if zcat fails to produce output use gunzip -c
                  if(@sso > 3000){ # if @sso is very big, I remove the useless contents
                      print "\n# open_sso_files: size of \@sso is for $file[$i] exceeds 3000 lines, ", scalar(@sso), " !!! \n";
                  }
                  push(@OUT, &read_sso_lines([@sso], $create_sso, $attach_range_in_names, $attach_range_in_names2,
                                  $new_format, $get_alignment) );
                  close SSO;
              }
         }
    }
    print "\n# \@OUT has ", scalar(@OUT), " elements \n" if $verbose;
	return(\@OUT); # @OUT has refs of hashes  (\%xxx, \%YYY, \%XXX,,,,)
}

#________________________________________________________________________
# Title     : write_fasta
# Usage     : many argments:  $seq_hash_reference  and $output_file_name
#             takes a hash which has got names keys and sequences values.
# Function  : writes multiple seqs. in fasta format (takes one or more seq.!!)
#             This needs hash which have 'name' 'actual sequence as value'
#
#             To print out each fasta seq into each single file, use write_fasta_seq_by_seq
#             This can rename seq names
#
# Example   : &write_fasta(\%in1, \$out_file_name, \%in2, \%in3,..., );
#             << The order of the hash and scalar ref. doesn't matter. >>
# Warning   : The default output file name is 'default_out.fa' if you do not
#             specify output file name.
#             OUTput file should have xxxxx.fa or xxxx.any_ext NOT just 'xxxxx'
# Keywords  : write_fasta_file, print_fasta_file, write fasta file, fasta_write
#             show_fasta
# Options   : v for STD out.
#             r for rename the sequences so that Clustalw would not complain with 10 char limit
#               so result wuld be:  0 ->ASDFASDF, 1->ASDFASFASF, 2->ADSFASDFA
# Returns   :
# Argument  :
# Version   : 2.3
#--------------------------------------------------------------------
sub write_fasta{
  #"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
  my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
  my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
  my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
  my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
  my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
  if($debug==1){print "\n\t\@hash=\"@hash\"
  \@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
  \@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  $| = 1;

  my($string, $string_leng, $na,$out_file_name_provided);
  my($output_file) ='default_out.fa'; ### when no output file name is given, this is used
  if(@file>0){
	$output_file = $file[0];
	$out_file_name_provided=1;
  }else{ $output_file='default_out.fa'; }

  for ($n=0 ; $n < @hash; $n ++){
	 my %hash=%{$hash[$n]};
	 my @keys=keys %hash;
	 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 # When only one seq is given, use the seq name as output file
	 #________________________________________________________________
	 if(@hash==1 and @keys==1 and @file < 1){
	    $output_file="$keys[0]\.fa";
	 }elsif(@file < 1){
	    $output_file="default_fa_$n\.fa";
	 }

	 open (FASTAS_WRITE,">$output_file");      # $string is the seq string.

	 for ($i=0; $i < @keys; $i++){
		$na= $keys[$i];
		$string = "\U$hash{$na}";
		$string=~s/[\n \.-]//g;	    # replaces all non-chars to null. '_' is used for stop codon
		if($char_opt=~/r/){  # rename the seqeunces with '0, 1, 2, 3," etc for  clustalw
		   $na=$i;
		}

		if($debug == 1){
			print ">$na\n";
			print FASTAS_WRITE ">$na\n";
	    }elsif($char_opt=~/v/){
		    print ">$na\n";
		    print FASTAS_WRITE ">$na\n";
		}else{
		    print  FASTAS_WRITE ">$na\n";
		}

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#  Main algorithm of writing in 60 char leng line
		#_____________________________________________________
		$string_leng=length($string);
		for($j=0; $j< $string_leng; $j+=60){
			if($debug == 1){
				printf "%.60s\n", substr($string,$j,60);
				printf FASTAS_WRITE "%.60s\n", substr($string,$j,60);
			}elsif($char_opt=~/v/i){
				printf "%.60s\n", substr($string,$j,60);
				printf FASTAS_WRITE "%.60s\n", substr($string,$j,60);
			}else{
			   printf FASTAS_WRITE "%.60s\n", substr($string,$j,60);
			}
		}
	 }
	 close FASTAS_WRITE;
  }
  if( $out_file_name_provided != 1){
	  print "\n\n# You didnt give out file name, $output_file  used\n";
  }
  if( -s $output_file ){
	 if($verbose=~/\S/){ ## if v option is given, mesg is omitted to prevent comments to a redirected output
	    print "\n# Sequences were written in  $output_file ";
	 }
  }else{
	 print "\n# The size of written outfile \"$output_file\" is 0, error \n\n"
  }
}
#________________________________________________________________________
# Title     : read_head_box
# Usage     : %entries = %{&read_head_box([\$file_to_read, \@BOXED ] )};
# Function  : Reads the introductory header box(the one you see on top of sub routines of
#             Jong's programs.). Make a hash(associative array) to put entries
#             and descriptions of the items. The hash values have new lines '\n' are
#             attached, so that later write_head_box just sorts Title to the top
#             and prints without much calculation.
#             This is similar to read_head_box, but
#             This has one long straight string as value(no \n inside)
#             There are two types of ending line one is Jong's #---------- ...
#             the other is Astrid's  #*************** ...
# Example   : Output is something like
#             ('Title', 'read_head_box', 'Tips', 'Use to parse doc', ...)
# Warning   :
# Keywords  : open_head_box, open_headbox, read_headbox
# Options   : 'b' for remove blank lines. This will remove all the entries
#             with no descriptions
# Returns   : A hash ref.
# Argument  : One or None. If you give an argu. it should be a ref. of an ARRAY
#              or a filename, or ref. of a filename.
#             If no arg is given, it reads SELF, ie. the program itself.
# Version   : 2.7
#--------------------------------------------------------------------
sub read_head_box{
  my($i, $c, $d, $j, $s, $z, @whole_file, $title_found, %Final_out,
	  $variable_string, $TITLE, $title, @keys, $end_found, $line, $entry,
	  $entry_match, $End_line_num, $remove_blank,  $title_entry_null,
	  $end_found, $Enclosed_entry, $Enclosed_var, $blank_counter,
	  $title_entry_exist, $entry_value, $temp_W, $Warning_part
	);

  if(ref($_[0]) eq 'ARRAY'){ ## When array is given
	  @whole_file = @{$_[0]};
  }elsif(-e ${$_[0]}){       ## When filename is given in a ref
	  open(FILE, "${$_[0]}");
	  @whole_file=(<FILE>);
  }elsif(-e $_[0]){          ## When filename is given
	  open(FILE, "$_[0]");
	  @whole_file=(<FILE>);
  }elsif( $_[0] eq 'b'){          ## When filename is given
	  $remove_blank = 1;
  }elsif( ${$_[0]} eq 'b'){          ## When filename is given
	  $remove_blank = 1;
  }else{
	  open(SELF, "$0");
	  @whole_file=(<SELF>);
  }
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  for($i=0; $i<@whole_file; $i++){
	 $whole_file[$i] =~ tr/\t/ {7}/;  ## This is quite important to some parsing!!!
	 #########################################
	 ##  The first and second line of box 1 ##
	 #########################################
	 if( ($whole_file[$i]=~/^#[_\*\~\-\=]{20,}$/)&&    ##  '#______' is discarded
		 ($whole_file[$i+1]=~/ *\# {0,3}([TitlNam]+e) {0,8}: {1,10}([\w\.:]*) *(Copyright.*)/i) ){
		 $TITLE = $1;      $title = "$2\n";   $Final_out{'Warning'}.="$3\n";
		 $entry_match=$TITLE; ## The very first $entry_match is set to 'Title' to prevent null entry
		 if($TITLE =~ /^Title|Name$/i){   #
			  if( ($title=~/^\s+$/)||( $title eq "\n") ){
				  $title_entry_null =1;  $title = '';  }    }
		 $Final_out{$TITLE}=$title;
		 $title_found ++ ;   $i++;  ## << this is essential to prevent reading the same line again.
		 last if $title_found > 1;    }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ## The first and second line of box 2, #__________ or #**************
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&($whole_file[$i]=~/^#[_\*]{20,}$/)&&
		 ($whole_file[$i+1]=~/^# {1,3}(\w{1,6}\s{0,2}\w+) {0,7}: {1,5}(.*) */i) ){
		 $title_found ++ ;        $i++;
		 $entry_match=$1;       $entry_value=$2;
		 $entry_match =~ s#^\S#(($tmp = $&) =~ tr/[a-z]/[A-Z]/,$tmp)#e;  ## Capitalize words
		 $Final_out{$entry_match}.= "$entry_value\n";
		 last if $title_found > 1;  next;   }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ##  'Enclosed' : section. After this, everything is read without discrimination ##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($Enclosed_entry == 1)&&($whole_file[$i] =~/^#{1,} {1,}(.*)$/) ){
		 $Final_out{$Enclosed_var}.= "$1\n";    }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ##  With proper entry 1 : for  'eg)'
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&($title_found==1)&&
		 ($whole_file[$i]=~ /^# {1,12}(eg ?\)) {0,8}(.*)/i)){
		 $entry_match='Example';
		 $Final_out{$entry_match}.= "$2\n";
	 }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ##  With PROPER entry 2 : descriptins like. 'Ussage : ssssssxxjkk  kj'
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&($title_found==1)&&
		 ($whole_file[$i]=~ /^# {1,2}(\w{1,4}\s{0,2}\w{1,7}) {0,8}[:\)] {0,6}(.*) */i)){
		 $entry_match=$1;       $entry_value=$2;
		 $entry_match =~ s#^\S#(($tmp = $&) =~ tr/[a-z]/[A-Z]/,$tmp)#e;
		 $Final_out{$entry_match}.= "$entry_value\n";
		 if($entry_match=~/^(Enclosed?)$/i){
			  $Enclosed_entry = 1;  $Enclosed_var=$1;        }    }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 #  With proper entry 3 : descriptins like. 'Ussage :', But blank description ##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&($title_found==1)&&
		 ($whole_file[$i]=~ /^# {1,2}(\w{1,4}\s{0,2}\w{1,7}) {0,8}[:\)]( {0,})$/i)){
		 $entry_match=$1;       $entry_value=$2;
		 $entry_match =~ s#^\S#(($tmp = $&) =~ tr/[a-z]/[A-Z]/,$tmp)#e;
		 $Final_out{$entry_match}.= " $entry_value\n";
		 if($entry_match=~/^(Enclosed?)$/i){
			  $Enclosed_entry = 1;  $Enclosed_var=$1;      }    }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 #  $option variable matching                ##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1) && ($title_found==1) &&
		 ($whole_file[$i]=~ /^# {1,15}([\$\@]\w+ +[\w=\>]+ +\S+ \w+ \S+ *.*)/ )){
		 $Final_out{$entry_match} .= "$1\n";  }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ###  all space line matching                 ##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&  ##<----- If blank line is matched. Take the line
		 ($title_found==1)&&($whole_file[$i]=~/^# {0,}$/) ){
		 $blank_counter++;
		 if($blank_counter > 2){ $blank_counter--; }
		 else{ $Final_out{$entry_match}.= " \n";  }     }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ###  Anything after 3 space to 12 positions  ##
	 ###  To match 'examples' etc. INC. ':'       ##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&
		 ($title_found==1)&&($whole_file[$i]=~/^#( {2,12})(.+)/) ){
		 $Final_out{$entry_match}.= "$1$2\n"; $blank_counter=0; }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ###  Anything after 1 space to 11 positions  ##
	 ###  To match 'examples' etc. EXC. ':'       ##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&
		 ($title_found==1)&&($whole_file[$i]=~/^# {1,12}([^:.]+)/) ){
		 $Final_out{$entry_match}.= "$1\n"; $blank_counter=0;}

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ###-------End of the read_box reading--------##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($title_found==1)&&
		 ($whole_file[$i]=~ /^#[\~\=\*\-]{15,}/)){  ## to match '#-----..' or '#******..'(Astrid's)
		 $End_line_num = $i;       $end_found++;
		 last;      }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 #  <<<<  Check if there is option table >>>>  #
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( (/^#{10,} option table of this program   #{10,}/)&&($end_found >=1) &&($title_found==1)){
		 $option_tb_found++; ### This is a global var.
	 }
  } ## < End of for loop


  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  ### If title is not there at all     ####
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  @keys=keys %Final_out;
  for(@keys){
	  if(/^Title$/i){    ## No Entry of Title at all??
		  $TITLE =$&;
		  $title_entry_exist = 1;
		  if($Final_out{$_}=~/^ *$/){   ## if Title => Null or just space
			  $title_entry_null = 1;    }  }  }

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  ### When title entry is not there    ####
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if( $title_entry_exist != 1){
		for($s=$End_line_num+1; $s < $End_line_num+20; $s++){
			if( $whole_file[$s] =~ /^sub {1,5}([\w\.]+) {0,6}\{/){
				$Final_out{'Title'} = "$1\n";   last;       }
			elsif( $whole_file[$s] =~/^#________________________+/){
				#######################################
				## Uses running file name as titile  ##
				#######################################
				$Final_out{'Title'} = "$0";     last;
			}else{
				#######################################
				## Uses running file name as titile  ##
				#######################################
				$Final_out{'Title'} = "$0";
			}
		}
  }
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  ### When title is blank              ####
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  elsif($title_entry_null ==1){  ## It looks for 'sub xxxxx{ ' line to get title
		### $End_line_num is the last line read.
		for($s = $End_line_num+1; $s < $End_line_num+20; $s++){
			if( $whole_file[$s] =~ /^sub {1,5}(\w+\.*\w*) {0,7}{/){
				$Final_out{$TITLE} = "$1\n";    last;     }
			elsif( $whole_file[$s] =~/^#________________________+/){
				#######################################
				## Uses running file name as titile  ##
				#######################################
				$Final_out{$TITLE} = "$0";     last;
			}else{
				#######################################
				## Uses running file name as titile  ##
				#######################################
				$Final_out{$TITLE} = "$0";
			}
		}
  }
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  ## Error handling, if no head box is found   ####
  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  if($title_found < 1){ print "\nFatal: No headbox found by read_head_box2 sub.\n";  }
  \%Final_out;
}               ##<<--- ENd of the sub read_head_box
#________________________________________________________________________
# Title     : open_fasta_files
# Usage     : %fasta_seq=%{&open_fasta_files($fasta_file, ['MJ0084'])};
#             if you put additional seq name as MJ0084 it will
#             fetch that sequence only in the database file.
#
#             %out=%{&open_fasta_files(@ARGV, \%index)};
#               while  %index has (seq indexpos seq2 indexpos2,,,)
#               In this case, the fasta file should have xxxx.fa format
#
# Function  : open fasta files and put sequences in a hash
#              If hash(es) is put which has sequence names and seek position
#              of the index file, it searches the input FASTA file to
#              fetch at that seek position. This is useful for Big fasta DBs
#             If the seq name has ranges like  XXXXXX_1-30, it will only
#              return 1-30 of XXXXXX sequence.
#
#             FASTA sequence file format is like this;
#
#             > 1st-seq
#             ABCDEFGHIJKLMOPABCDEFGHIJKLMOPABCDEFGHIJKLMOPABCDEFG
#             > 2nd.sequ
#             ABCDEFGHIJKLMOYYUIUUIUIYIKLMOPABCDEFGHIJKLMOPABCDEFG
#             >owl|P04439|1A03_HUMAN HLA CLASS I HISTOCOMPATIBILITY ANTIGEN, A-3 ALPHA CHAIN PRECURSOR....
#             MARGDQAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDT
#
#             This can also return the sizes of sequences rather than seqs.
#
#             This ignores any dup entrynames coming later.
#
# Example   : %out = %{&open_fasta_files(@ARGV)};
#             %out2=%{&open_fasta_files('seq.fa', \%index)};
#             %out3=%{&open_fasta_files('seq.fa', \%range)};
#             %seq=%{&open_fasta_files($PDB40_FASTA, \@seq_to_fetch)};
#
#             while @ARGV at prompt was: 'GMJ.pep MJ0084'
#
# Keywords  : open_fasta, open_fa_files, open_FASTA_files,
# Options   : Seq name to fetch the specified seq only.
#             as open_fasta_files.pl MY_SEQ_NAME Swissprot.fasta
#            -d  for giving back desc as well as the name. so it
#                gives  'HI0002 This is the description part'
#                as the key
#             If you put hash which is like('seq_name', ['20-30', '30-44',..])
#              it will produce hash which has got:
#              ( seq_name_20-30 'asdfasdfasdfasdfasd',
#                seq_name_30-44 'kljkljkjkjljkjljkll',
#                ....           .... )
#            -s for returning sequence size only
# Version   : 3.9
#--------------------------------------------------------------------
sub open_fasta_files{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

   my (%sequence, %HASH, @Keys, $seq_found1, $S_start, $S_end, $seq_found,
	   $present_seq, @seq_Names, %Sizes, $bare_seq_name, $fasta_seq_idx_file,
	   %seq_fragments);

   if(@file<1){
	  print "\n \@file has less than 1 elem. There is no fileinput for open_fasta_files\n";
	  exit
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  (1) When perl file INDEX pos info is given in hash, this speeds up
   #_________________________________________________________________________________
   for($d=0; $d < @hash; $d++){
	   my ($sequence, $NAME, $range_start, $range_leng);
	   %HASH=%{$hash[$d]};
	   my @Keys=keys %HASH;  ## <<< NOTE it is @Keys, not @keys
	   for($f=0; $f< @file; $f++){
		  #====== It must be xxxx.fa format =======
		  unless($file[$f]=~/\S\.fa[sta]?$/){
			  print "\n# open_fasta_files: \$file\[\$f\] does not have fasta extension, skipping"; next; }
		  open(FASTA, $file[$f]);
		  F0: for($e=0; $e< @Keys; $e++){
			 my $sequence;
			 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			 # When seq name has range attachment, it handles
			 #________________________________________________
			 if($Keys[$e]=~/^(\S+)_(\d+)\-(\d+)/){
				 $NAME=$1;
				 $range_start=$2-1;    ## to fit in substr function
				 $range_leng =$3-$2+1; ## to fit in substr
			 }else{
			     $NAME=$Keys[$e];
			 }
			 if($HASH{$Keys[$e]}=~/^(\d+)$/){
				 splice(@hash, $d, 1);
				 $d--;
				 splice(@file, $f, 1);
				 $f--;
				 seek(FASTA, $1-220, 0);  # -220 is necessary
				 while(<FASTA>){
					 if( /^\> *$NAME/  or
						 /^\> *owl\|\S+\|$NAME/){  # to handle ">owl|P04439|1A03_HUMAN HLA CLASS I HISTOCOMPATIBILITY
					        $seq_found1=1;
					 }elsif(/^(\w+)$/ and $seq_found1==1){	 $sequence .=$1;
					 }elsif(/^\> *\S+/ and $seq_found1==1){
						  #======= When range is defined, take only the ranged part==================
						  if($range_start =~/\d+/){
							  $sequence{$Keys[$e]}=substr($sequence, $range_start, $range_leng);
						  }else{	 $sequence{$Keys[$e]}=$sequence; }
						  $range_start='';
						  $sequence='';
						  $seq_found1=0; next F0;
					 }
				 }
			  }
		  }
	  }
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
   # (2) opening FASTA files (NORmal, when no perl index pos number is given)
   #_______________________________________________________________________
   for($i=0; $i< @file; $i++){
	   unless(-s $file[$i]){ next; } ## this is essential as handle_arguments has a problem
	   my($entry_found, $name, $matched);
	   my($input_file) = ${$file[$i]} || $file[$i];

	   if($debug eq 1){ print "\n open_fasta_files: Inputfile is $input_file\n" };
	   unless (-e $input_file){
			print chr(7);
			print "\n\n\t This is sub open_fas_files in $0  \n\n";
			print "\t Fatal: The input file $input_file is not in the directory \n";
	   }
	   open(FILE_1,"$input_file");
	   if(@hash >=1){  ## if seq names are given in hash
		   for($h=0; $h< @hash; $h++){
			  @string=(@string, keys %{$hash[$h]});
		   }
	   }
	   @string=sort @string;
	   $num_of_seq_to_fetch=@string;
	   if(@string > 0){
		   print "\n# open_fasta_files(normal fasta fetch): \$num_of_seq_to_fetch is $num_of_seq_to_fetch\n";
	   }

	   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	   #  (2.1) when seq to fetch is given by \@sequences  format
	   #_______________________________________________________________________
	   if( @_ > 1  and  @string > 0 ){
		   print "\n#  open_fasta_files is fetching sequences from \$input_file= $input_file\n";
		   %sequence=%{&fetch_sequence_from_db($input_file, \@string)};
		   print "\n# $fasta_seq_idx_file file is made by open_fasta_files(fetch_sequence_from_db), you may remove it\n";
	   }
	   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	   #  (2.2) When seq names NOT given, fetches all (THE DEFAULT)
	   #____________________________________________________________
	   else{
		 while(<FILE_1>){                # file1 needs to be xxxx.fasta for the moment, automatic later
			if(/^> *gi\|\d+\|\S+\|(\S+)\|.*/){  ## for >gi|1669546|dbj|D84107|D84107 Human mRNA for Werner syndrome-1/type 1, complete cds
				 if($char_opt=~/[\-]?d/i){  # To add the description
					 $name=$_;  # entire line becomes the name of the seque.
				 }else{
					 if( $sequence{$1} ){
						 #------- To avoid identical entry reading repeatedly -----
						 print "\n# I am open_fasta_files: $1 seems to be the same as previous entry, ERROR??\n";
						 $entry_found=0;
					 }else{      $name=$1;   $entry_found=1;     }
				 }
			}elsif(/^> *owl\|\S+\|(\S+)/){  ## for ">owl|P04439|1A03_HUMAN HLA CLASS I HISTOCOMPATIB
				 if($char_opt=~/[\-]?d/i){  # To add the description
					 $name=$_;  # entire line becomes the name of the seque.
				 }else{
					 if( $sequence{$1} ){
						 #------- To avoid identical entry reading repeatedly -----
						 print "\n# I am open_fasta_files: $1 seems to be the same as previous entry, ERROR??\n";
						 $entry_found=0;
					 }else{      $name=$1;   $entry_found=1;     }
				 }
			}elsif(/^> {0,5}([\w\-\.]+) *.*$/){
				 if($char_opt=~/[\-]?d/i){   $name=$_;  # To add the description
				 }else{
					 if( $sequence{$1} ){ # check if the entry already exists
						print "\n# $1 seems to be the same as previous entry, ERROR??\n";
						$entry_found=0;
					 }else{     $name=$1;   $entry_found=1;      }
				 }
			}elsif(/^([\w\.\- ]+)$/ and $entry_found == 1){
				 $matched=$1;    $matched=~s/ //g;
				 $sequence{$name}.= $matched if defined($name);
			}elsif(/^$/){  next;
			}else{  $entry_found=0;  } ## this is when rubbish is matched
		 }# end of while
	   }
	   close FILE_1;
   }


   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~~~~~~~~~~~~~~~~~~~~~`
   # (3) When ranges information is given(via \@range), seq in those ranges are returned
   #______________________________________________________________________________________
   if(defined(@range)){
	   %seq_fragments=%{&get_seq_fragments(\%sequence, \@range)};
	   return(\%seq_fragments);
   }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
   # (4) When only size is asked with -s option
   #_____________________________________________________________________________
   elsif($char_opt=~/s/){ # when SIZE(length of seq) return only option is set
	   @seq_Names=keys %sequence;
	   for($i=0; $i<@seq_Names; $i++){
		  $Sizes{$seq_Names[$i]}=length($sequence{$seq_Names[$i]});
	   }
	   return(\%Sizes);
   }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
   # (5) when hash which has range info is given(@range should not be defined)
   #_____________________________________________________________________________
   elsif(@hash >=1){
	   for($h=0; $h< @hash; $h++){
		   my %hash=%{$hash[$h]};
		   my @Keys=keys %hash;
		   for($k=0; $k<@Keys; $k++){
			   if(defined($hash{$Keys[$k]})){
				  ($S_start, $S_end)=$hash{$Keys[$k]}=~/(\d+)\-(\d+)/;
				  $sequence{$Keys[$k]}=substr($sequence{$Keys[$k]}, ($S_start-1), ($S_end-$S_start));
			   }
		   }
	   }
	   return(\%sequence);
   }else{
	   return(\%sequence);
   }
}
#______________________________________________________________________________
# Title     : make_seq_index_file
# Usage     : @idx_files_made=@{&make_seq_index_file(\@file)};
# Function  : creates xxxx.fa.idx file and makes a link to pwd. If @file contains
#              names with .idx extension already, it will not put another idx
#              index to it.
# Example   :
# Keywords  : make_fasta_seq_index_file, create_seq_index_file, make_idx_file,
#             create_idx_file, create_seq_idx_file, make_index_file, create_index_file
#             make_sequence_index_file, create_sequene_index_file
# Options   :
# Version   : 1.3
#----------------------------------------------------------------------------
sub make_seq_index_file{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	my(@index_files_made, $fasta_db_input, $fasta_db_idx, %index);
	print "\n# make_seq_index_file : input \@file was @file\n";

	for($i=0; $i< @file; $i++){
		$fasta_db_input=$file[$i];
		if($fasta_db_input=~/\S+\.idx$/){
		    $fasta_db_idx=$fasta_db_input;
		}else{
			$fasta_db_idx="$fasta_db_input.idx";
		}

		open(FASTA_DB, "$fasta_db_input");
		open(FASTAIDX, ">$fasta_db_idx");

 	    print FASTAIDX "# fasta_index for $fasta_db_input\n";
		while(<FASTA_DB>){
			if(/^\> {0,4}(\S+) */){
				$index{$1}=tell(FASTA_DB);
				print FASTAIDX "\n$1 $index{$1}";
			}
		}
		close(FASTA_DB, FASTAIDX);
		if(-s $fasta_db_idx){
			print "\n# The size of $fasta_db_idx is more than 0, looks O.K. \n";
			push(@index_files_made, $fasta_db_idx);
			system("ln -s $fasta_db_idx .");
		}else{
		    print "\n# The size of $fasta_db_idx is less than 0, ERROR??\n";
		}
	}
	if(@file < 2){
	   return( \$fasta_db_idx );
	}else{
	   return(\@index_files_made);
	}
}
#_____________________________________________________________________________
# Title     : fetch_sequence_from_db
# Usage     : %sequence=%{&fetch_sequence_from_db($input_file, \@string)};
# Function  : accept seq names (with or without ranges like _10-111 )
#              and produces hash ref.
#             As an option, you can write(xxxx.fa) the sequences in pwd
#              with the file names with sequence names.
#             The default database used is FASTA format OWL database.
#              You can change this by S (for Swissprot either fasta
#              or full format), P for PDB fasta format data.
#             If you give the path name of DB, it will look for the
#              DB given.
#
#             This automatically checks sequence family number as
#               in >d1bpi___7.6.1
#               and attaches the number in final %sequence output
#
# Example   : %seq=%{&fetch_sequence_from_db(\@input, seq.fa, seq.fa.idx)};
#              while @input=qw( 11S3_HELAN_11-31 A1AB_CANFA A1AT_PIG )
# Keywords  : fetch_seq_from_db, fetch_sequence_from_database
# Options   : _  or #  for debugging.
#     w       for write fasta file
#     d=p100  for PDB100 fasta database from ENV
#     d=p40   for PDB40  fasta database from ENV
#     d=p     for PDB database (usually p100) from ENV
#     d=s     for Swissprot database from ENV
#     d=o     for OWL database from ENV
#     i=      for index filename. If not specified, this looks for it in the same dir as fast     ò
#          ¨√@ê
# Returns   : ref of hash
# Argument  : gets names of sequences
#             eg) \@array, \%hash, \$seq, while @array=(seq1, seq2), $seq='seq1 seq1'
#                                               %hash=(seq1, xxxx, seq2, yyyy);
# Version   : 2.9
#------------------------------------------------------------------------------
sub fetch_sequence_from_db{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

	my(@DATABASE, @INDEX_FILE, %sequence, %seq_with_index, @input_seq_names,
	   %long_index, @Keys, $R_start, $NAME, $R_leng, $found_seq_count,
	   $seq_found1, $sequence, @keys, $index_file);

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# getting input seq names from all sources
	#________________________________________________________
	for(0..@hash){
	    push(@input_seq_names, keys %{$hash[$_]} );
	}
	for(0..@raw_string){
		push(@input_seq_names, split(/ +/, $raw_string[$_]) );
	}
	print "\n# (1) fetch_sequence_from_db: \@raw_string has: ", scalar(@raw_string), " elements";
	print "\n# (2) fetch_sequence_from_db: No. of seq to fetch is:",scalar(@input_seq_names);
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Choose the DBs and INDEX for fetching sequences. All files input must be DATABAE or INDEXfile
	#___________________________________________
	if(@file > 0){
		for($i=0; $i< @file; $i++){
		   if(-T $file[$i] and $file[$i]=~/\.fa[sta]?$/){      push(@DATABASE, $file[$i]);   next}
		   elsif((-T $file[$i]) and ($file[$i]=~/\.seq$/)){    push(@DATABASE, $file[$i]);   next}
		   elsif((-T $file[$i]) and ($file[$i]=~/\.dat$/)){    push(@DATABASE, $file[$i]);   next}
		   elsif(-T $file[$i] and $file[$i]=~/\.idx$/){        push(@INDEX_FILE, $file[$i]); next }
		   if($file[$i] !~/\.idx/ and -T "$file[$i]\.idx"){    push(@INDEX_FILE, "$file[$i]\.idx"); }
		   else{
			  print "\n#  WARN:  fetch_sequence_from_db:
			  You put a file-name-like which is not a fasta DB. Error. I am removing $file[$i]";
			  splice(@file, $i, 1);
			  $i--;
		   }
		}
	}

	if($vars{'d'}=~/^p[100]*$/){
	   if( -T  $ENV{'PDB_FASTA'} ){             push(@DATABASE,   $ENV{'PDB_FASTA'} );     }
	   elsif(  -T $ENV{'PDB_SEQ_FASTA'} ){      push(@DATABASE,   $ENV{'PDB_SEQ_FASTA'}  ); }
	   elsif(  -T $ENV{'PDB100_FASTA'} ){       push(@DATABASE,   $ENV{'PDB100_FASTA'} ); }
	   if(  -T $ENV{'PDB_FASTA_INDEX'} ){       push(@INDEX_FILE, $ENV{'PDB_FASTA_INDEX'} ); }
	}elsif( $vars{'d'}=~/^p\d+d$/ ){
	   if(  -T $ENV{'PDB100D_FASTA'} ){         push(@DATABASE,   $ENV{'PDB100D_FASTA'});     }
	   elsif(  -T $ENV{'PDBD100_FASTA'}  ){     push(@DATABASE,   $ENV{'PDBD100_FASTA'}); }
	   elsif(  -T $ENV{'PDB100D_SEQ_FASTA'}  ){ push(@DATABASE,   $ENV{'PDB100D_SEQ_FASTA'}); }
	   elsif(  -T $ENV{'PDBD100_SEQ_FASTA'}  ){ push(@DATABASE,   $ENV{'PDBD100_SEQ_FASTA'}); }
	   if(  -T $ENV{'PDB100D_SEQ_FASTA_INDEX'} ){    push(@INDEX_FILE, $ENV{'PDB100D_SEQ_FASTA_INDEX'}); }
	   elsif(  -T $ENV{'PDBD100_SEQ_FASTA_INDEX'} ){    push(@INDEX_FILE, $ENV{'PDBD100_SEQ_FASTA_INDEX'}); }
	}elsif( $vars{'d'}=~/^p40/ ){
	   if(  -T $ENV{'PDB40_FASTA'} ){          push(@DATABASE,   $ENV{'PDB40_FASTA'});     }
	   elsif(  -T $ENV{'PDB40_SEQ_FASTA'}  ){  push(@DATABASE,   $ENV{'PDB40_SEQ_FASTA'}); }
	   if(  -T $ENV{'PDB40_FASTA_INDEX'} ){    push(@INDEX_FILE, $ENV{'PDB40_FASTA_INDEX'}); }
	}elsif( $vars{'d'}=~/^p90/ ){
	   if(  -T $ENV{'PDB90_FASTA'}  ){         push(@DATABASE,   $ENV{'PDB90_FASTA'}    ); }
	   elsif(  -T $ENV{'PDB90_SEQ_FASTA'} ){   push(@DATABASE,   $ENV{'PDB90_SEQ_FASTA'}); }
	   if(  -T $ENV{'PDB90_FASTA_INDEX'} ){    push(@INDEX_FILE, $ENV{'PDB90_FASTA_INDEX'}); }
	}
	if( $vars{'d'}=~/s/){
	   if(  -T $ENV{'SWISS_FASTA'} ){          push(@DATABASE,   $ENV{'SWISS_FASTA'});     }
	   elsif(  -T $ENV{'SWISS_SEQ_FASTA'} ){   push(@DATABASE,   $ENV{'SWISS_SEQ_FASTA'}); }
	   elsif(  -T $ENV{"SWISS_DIR\/seq.fa"} ){ push(@DATABASE,   $ENV{"SWISS_DIR\/seq.fa"}); }
	   if(  -T $ENV{'SWISS_FASTA_INDEX'} ){    push(@INDEX_FILE, $ENV{'SWISS_FASTA_INDEX'}); }
	   elsif(  -T $ENV{'SWINDEX'} ){           push(@INDEX_FILE, $ENV{'SWINDEX'}); }
	}
	if( $vars{'d'}=~/^o/){
		if(  -T $ENV{'OWL_FASTA'} ){            push(@DATABASE,   $ENV{'OWL_FASTA'});     }
		elsif(  -T $ENV{'OWL_SEQ_FASTA'} ){     push(@DATABASE,   $ENV{'OWL_SEQ_FASTA'}); }
		elsif(  -T $ENV{"OWL_DIR\/owl.fa"} ){   push(@DATABASE,   $ENV{"OWL_DIR\/owl.fa"}); }
		if(  -T $ENV{'OWL_FASTA_INDEX'} ){      push(@INDEX_FILE, $ENV{'OWL_FASTA_INDEX'}); }
		print "\n# Fetching sequences from OWL\n";
	}
	if( $vars{'d'}=~/^\S+\.\S+$/ and -T $vars{'d'} ){ push(@DATABASE, $vars{'d'} );     }
	if( $vars{'i'}=~/\S+\.\S+$/ and -T $vars{'i'} ){ push(@INDEX_FILE, $vars{'i'} );   }
	if(@INDEX_FILE > 0 and @DATABASE > 0){
		if( ${&if_file_older_than_x_days("$DATABASE[0]\.idx", 5)} > 0 ){
			$index_file=${&make_seq_index_file(\@DATABASE)};
	        push(@INDEX_FILE, $index_file);
		}elsif(-s "$DATABASE[0]\.idx"){
			push(@INDEX_FILE, "$DATABASE[0]\.idx");
		}else{
			print "\n# fetch_sequence_from_db: Some weird error in pushing \$index_file to \@INDEX_FILE\n"; exit;
		}
	}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#  Final check for ALL the inputs
	#___________________________________________________
	if( @DATABASE  < 1){ print "\n# fetch_sequence_from_db: DATABASE file no found. Error\n"; exit     }
	if( @INDEX_FILE < 1){
		print "\n# fetch_sequence_from_db: \@INDEX_FILE has less than 1 elem. Error\n";
		push(@INDEX_FILE, ${&make_seq_index_file(@DATABASE)});
		print "     fetch_sequence_from_db called make_seq_index_file to make @INDEX_FILE\n";
	}
 	if($debug==1){
		print "\n# DATABASE used     : @DATABASE";
		print "\n# INDEX_FILE used   : @INDEX_FILE";
		print "\n# input_seq_names   : @input_seq_names";
	}

	#--------------------------------------------------------------------------


  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##  Now I have @DATABASE, @INDEX_FILE, @input_seq_names
  ##_______________________________________________________________

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	#  Reading in index file to get 'seq' 'seek pos' to make %seq_with_index
	#__________________________________________________________________________
	print "\n#  fetch_sequence_from_db: \@INDEX_FILE @INDEX_FILE, \@DATABASE :@DATABASE\n";
	for($i=0; $i< @INDEX_FILE; $i++){
	   open(INDEX, "$INDEX_FILE[$i]");
	   while(<INDEX>){
		  if(/(\S+) +(\S+)/){
			  $long_index{$1}=$2;
		  }
	   }
	   for($j =0; $j < @input_seq_names; $j++){

			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~``
			#  If DATABASE has sequence names with ranges already index the seq with ranges
			#____________________________________________________________________________________
			if($input_seq_names[$j]=~/(\S+\_\d+\-\d+)/ and $long_index{$1}){

			    $seq_with_index{$1}=$long_index{$1};

			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~``
			#  If DATABASE has sequence names without ranges index the seq without ranges
			#____________________________________________________________________________________
			}elsif($input_seq_names[$j]=~/(\S+)\_\d+\-\d+/ and $long_index{$1}){

				$seq_with_index{$input_seq_names[$j]}=$long_index{$1}; # !!!! <--- This line is critical

			}elsif($input_seq_names[$j]=~/(\S+)\_\d+\-\d+/ and $long_index{"$1\_"}){ # to handle Tim's new pdb100.fa files

			    $seq_with_index{$input_seq_names[$j]}=$long_index{"$1\_"};
			    print "\n# Warning: $1 (from $input_seq_names[$j]) matched with $1\_ in $INDEX_FILE[$i],
					  I hope this is correct!!\n";
			}
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~``
			#  If input_seq_name has SCOP superfamily numbers
			#____________________________________________________________________________________
			elsif($input_seq_names[$j]=~/^(\S+)\_(\d+\.\d+\.\d+)[\.\d+\.\d+]*/ and $long_index{$1}){

				$seq_with_index{"$1\_$2"}=$long_index{$1}; # !!!! <--- This line is critical

			}elsif($input_seq_names[$j]=~/\S/ and $long_index{$input_seq_names[$j]}){
				$seq_with_index{$input_seq_names[$j]}=$long_index{$input_seq_names[$j]}
			}else{
				print "\n#  $input_seq_names[$j](with, without range) have NO corresponding index in $INDEX_FILE[$i], ERR";
			}
	   }
	   close INDEX;
	   if ( scalar(keys %seq_with_index) < 1){
		    print "\n# fetch_sequence_from_db: \%seq_with_index is too small, ERROR?\n";
	   }
	}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
	#  Fetching sequences from DATABASE
	#_______________________________________________________________
	print "\n# fetch_sequence_from_db: Fetching seqs from @DATABASE with  @INDEX_FILE ";
	my @Keys= keys %seq_with_index;        ## <<< NOTE it is @Keys, not @keys
	print "\n# (3) fetch_sequence_from_db: No. of seq indexed is:", scalar(@Keys);

	for($f=0; $f< @DATABASE; $f++){
	   open(FASTA, $DATABASE[$f]);
	   F0: for($e=0; $e< @Keys; $e++){
		  my ($seq_found1, $super_fam_class, $NAME, $R_leng, $R_start, $sequence);
		  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		  # When seq name has range attachment, it handles
		  #________________________________________________
		  if($Keys[$e]=~/(\S+)_(\d+)\-(\d+)$/){
			  $NAME=$1;
			  $R_start=$2-1;      ## to fit in substr function
			  $R_leng =$3-$2+1; ## to fit in substr
			  print "\n# (4) fetch_sequence_from_db: Sequences have ranges only (not superfamily numb.) \n";
		  }
		  elsif($Keys[$e]=~/(\S+)_(\d+)\-(\d+)\_(\d+\.\d+\.\d+)[\.\d+\.\d+]*/){
			  $NAME=$1;
			  $R_start=$2-1;      ## to fit in substr function
			  $R_leng =$3-$2+1; ## to fit in substr
			  $super_fam_class=$4;
			  print "\n# (4) fetch_sequence_from_db: Sequences have ranges and superfamily numb.\n";
		  }
		  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		  # When superfamily (scop) number is attached
		  #___________________________________________________
		  elsif($Keys[$e]=~/(\S+)\_(\d+\.\d+\.\d+)[\.\d+\.\d+]*/){
			  $NAME=$1;
			  $super_fam_class=$2;
			  print "\n# (4) fetch_sequence_from_db: Sequences have SCOP superfamily numbers only \n";
		  }elsif($Keys[$e]=~/^ *(\S+)[\,]*$/){
			  print "\n# (4) fetch_sequence_from_db: Sequences DON't have ranges or SCOP superfam numb.\n";
			  $NAME=$1;
		  }

		  if($seq_with_index{$NAME}=~/(\d+)/        # It is importnt having $seq_with_index{$Keys[$e]}
			   or $seq_with_index{$Keys[$e]}=~/(\d+)/
			   or $seq_with_index{"$NAME\,"}=~/(\d+)/    # this is for overcoming '>xxxx,'  entry(the comma)
			   or $seq_with_index{"$NAME\_"}=~/(\d+)/    # to handle Tim's  >c1eru_ 3.30.1.1.4
			   ){
			   my $finding_position= $1-300;
			   if( $finding_position >= 0 ){   seek(FASTA, $1-300, 0);  # -300 is necessary
			   }elsif($finding_position < 0){  seek(FASTA, 0, 0); }      ## This is essential !!!
			   while(<FASTA>){
				  if(!$seq_found1){
					  if(/\> *$NAME[\,_]? +\d+\./){
						  $seq_found1=1;
					  }
				  }else{
					  if(/ *(\w+) *$/ ){ $sequence .=$1;  ## you should use $1 to avoid including NEW line
						  unless(eof FASTA){ next   ## This is critically important to prevent error.
						  }else{ goto PUT_SEQ }     ## If the last seq has only one single line seq string, it could be a problem
					  }elsif( (/ *\> *\S+/)  or (eof FASTA) ){
						  #======= When range is defined ==================
						  PUT_SEQ:
						  if($R_start =~/\d+/){
							  $sequence{$Keys[$e]}=substr($sequence, $R_start, $R_leng); next; #
						  }
						  #======= To handle superfamily information given ==========
						  if($super_fam_class){
							  $sequence{$Keys[$e]}=$sequence;
							  $acquired_seq_count++;
						  }
						  #======= When range is NOT defined ==================
						  else{
							  $sequence{$Keys[$e]}=$sequence;
						  }
						  ($R_start, $sequence, $seq_found1)='';  ## reset $R_start, $seq_found1,,
						  next F0;
					  }
				  }
			   }

		  }else{
			   print "\n# Error, the sequence pos for $NAME (from $Keys[$e]) in DB doesnt exist in xxxx.idx file?\n";
		  }
	   }
	   close FASTA;
	}
	#print "\n# (6) fetch_sequence_from_db: counted fetched seqs: $found_seq_count, $acquired_seq_count";
	#print "\n# (7) fetch_sequence_from_db: Fetching seq has finished \n";

	return(\%sequence);
}
#________________________________________________________________________
# Title     : reverse_sequences
# Usage     : %out = %{&rev_sequence_one_hash(\%input_seq_hash, \%hash2,...)};
# Function  : gets ref. of strings, reverses the elems.
# Example   :
# Warning   :
# Keywords  : reverse_sequence, reverse_sequence_hash, rev_sequence_hash
# Options   :
# Returns   : one or more hash references.
# Argument  : hash, eg(1, 'skdfj', 2, 'kdfjkdj', 3, 'kdfjk');
#             Input example:
#             ..
#             >HI0256
#             FLSANVLPIAPIINGGRTAVDNITQSVSDKPFVKDIGTKIKEAIALSKYSTQPQYISTTN
#             >HI0094
#             DILRTFVKMETGLKFPKKFKLKANLALFMNRRNKRPDTIMTAVADAGQKISEAKLNTTAK
#             ..
#
#             Output example: (Reversed :-)
#             ..
#             >HI0256_rv   <<-- note the added extension
#             ALDJFLKAJFJALSDJFLAJSLFJAKLSDFJLASJDFLAJSLDFJASJDFLJSDFJSDLJ
#             >HI0094_rv
#             LASJDFLKAJFJALSDJFLKSDJLFAJLKDJFLASJDFLKDFJKDJFKDJFKDJFKJDLJ
#             ..
#
# Version   : 1.2
#--------------------------------------------------------------------
sub reverse_sequences{
	my(%rev_hash, @rev_hash_refs, $name, $name_with_ext, $i);
	for($i=0; $i < @_; $i++){
	    my %in_hash = %{$_[$i]};
		my @keys    = keys %in_hash;
		for $name (@keys ){
		    $name_with_ext = "$name\_rv";
			$rev_hash{$name_with_ext} = reverse($in_hash{$name});
		}
		push(@rev_hash_refs, \%rev_hash);
	}
	if(@rev_hash_refs ==1){ return($rev_hash_refs[0]);}
	else{ return(@rev_hash_refs);}
}
#_________________________________________________________________________________
# Title     : read_sso_lines
# Usage     : &read_sso_lines([@sso], $create_sso, $attach_range_in_names, $attach_range_in_names2,
#                  $new_format, $get_alignment) );
# Function  : Main subroutine for open_sso_files.
# Example   :
# Keywords  : read_sso_lines_in_array
# Options   : a c r r2 n
# Version   : 1.0
#----------------------------------------------------------------------------
sub read_sso_lines{
	  my ($i, $upper_expect_limit, $lower_expect_limit)=(50,0); ##<<--- DEFAULT
	  my (@out_refs, $parseable, @SSO, $create_sso, $i, $j, $k, $attach_range_in_names);

	  for($i=0; $i< @_; $i++){
		  if($_[$i]=~/u=(\S+)/){    $upper_expect_limit=$1 }
		  elsif(ref($_[$i]) eq 'ARRAY'){ @SSO=@{$_[$i]};   }
		  elsif($_[$i]=~/l=(\S+)/){ $lower_expect_limit=$1 }
          elsif($_[$i]=~/^c *$/){     $create_sso = 'c' ; print "\n# read_sso_lines: \$create_sso is set";}
		  elsif($_[$i]=~/^a *$/){     $get_alignment='a'; }
		  elsif($_[$i]=~/^r *$/){   $attach_range_in_names='r' }
		  elsif($_[$i]=~/^r2 *$/){   $attach_range_in_names2='r2' }
		  elsif($_[$i]=~/^n *$/){   $new_format='n' }
	  }
	  print "\n# \$attach_range_in_names2 is $attach_range_in_names2\n" if $attach_range_in_names2;

	  #~~~~~~ Checking if sso is a parseable form or not~~~~~~~~~~~~~
	  TEMP:for($k=0; $k < @SSO; $k++){
		  if($SSO[$k] =~ /\>\>\>/  or $SSO[$k] =~ /^ *\; \S+\:/ ){
			  $parseable++;  if($parseable >= 10){  last TEMP;     }
		  }elsif($SSO[$k]=~/^  +\:+/){ $parseable--;
		  }elsif($SSO[$k] =~ /^ +1\>\>\>(\S+)/){ $pvm_version_fasta_out=1; $parseable +=10; $original_target_seq=$1; last TEMP;
		  }
	  }
	  if($parseable >= 10){
          @out_refs=@{&read_machine_readable_sso_lines(\@SSO, $get_alignment, $create_sso, $upper_expect_limit,
						   $new_format, $lower_expect_limit, $attach_range_in_names, $attach_range_in_names2)};
	  }else{
		  @out_refs=@{&read_machine_unreadable_sso_lines(\@SSO, $get_alignment, $create_sso, $upper_expect_limit,
						   $new_format, $lower_expect_limit, $attach_range_in_names, $attach_range_in_names2)};
	  }
	  return(@out_refs);
}
#___________________________________________________________________
# Title     : scramble_array
# Usage     : @in=@{&scramble_array(\@in)};
# Function  : shuffles the elements of array
# Example   :
# Keywords  : randomise_array, randomize_array, shuffle_array
# Options   :
# Version   : 1.4
#---------------------------------------------------------------
sub scramble_array{
	srand(time()|$$);  # or use srand(time^$$);
	my ($i, @scrambled, @out, @each_array);

	for($i =0; $i< @_; $i++){
	   my @each_array = @{$_[$i]};
	   while (@each_array) {
		   push @scrambled, splice @each_array, int(rand(@each_array)), 1;
	   }
	   push(@out, \@scrambled);
	}
	if(@out > 1){
	   return(@out);
	}else{
	   return($out[0]);
	}
}
#________________________________________________________________________
# Title     : write_fasta_seq_by_seq
# Usage     : &write_fasta_seq_by_seq(\%hash, [$extension], [\$output_filename]);
# Function  : accepts one hash of multiple sequences and writes many files
#             of single sequences by using the names as file names.
#             If $extension is provided, it writes an output as in
#             the below example (seq1_sc.fasta). If not, it just attach
#             'fa' to files.
#             This needs, hash of 'name', 'actual sequence as value'
# Example   : with >xxxx
#                  ASDFASDFASDFASDFASDFASDFASDF
#                  >yyyy
#                  ASDFASDFASDFASDFASDFASDFSDAFSD
#
#             You will get two files (xxxx.fa, yyyy.fa)
# Keywords  : write_each_fasta, write_single_fasta, write_fasta_single
#             single_fasta_write, write_fasta_files_seq_by_seq,
#             write_single_fasta_files,
# Options   : can specify extension name.
#             e  for checking fasta file exists or not and skipps if so
#             r for rename the sequences so that Clustalw would not complain with 10 char limit
#               so result wuld be:  0 ->ASDFASDF, 1->ASDFASFASF, 2->ADSFASDFA
# Returns   : nothing. default OUTPUT file name is '$key.fa' !!
# Version   : 1.9
#--------------------------------------------------------------------
sub write_fasta_seq_by_seq{
	 my ($i, $exists_opt, $rename_seq_opt, $out_file_name_given);
	 for($i=0; $i< @_; $i++){
		if($_[$i]=~/e$/){
		   $exists_opt=1;
		   splice(@_, $i, 1);
		   $i--;
		}elsif($_[$i]=~/r$/){
		   $rename_seq_opt='r';
		   splice(@_, $i, 1);
		   $i--;
		}elsif( $_[$i] =~/\.fa/ or -e $_[$i] ){
		   $out_file_name_given=1;
		   $out_file_name = $_[$i];
		   splice(@_, $i, 1);
		   $i--;
		}elsif( ref ($_[$i]) eq 'SCALAR'){
		   if( ${$_[$i]} =~/\.fa/ or -e ${$_[$i]} ){
		      $out_file_name_given=1;
		      $output_file=${$_[$i]};
		      splice(@_, $i, 1);
		      $i--;
		   }
		}
	 }
	 my(%temp_hash, $key, $output_file);
	 my(%input)     =%{$_[0]};
	 my($extension) =${$_[1]} || $_[1];
	 for $key (keys %input){
		my %temp_hash=();
		$temp_hash{$key}=$input{$key};
		if (($key=~ /\_$extension$/)||($#_ == 0)){
			$output_file = "$key\.fa";
		}else{
			$output_file = "$key\_$extension\.fa";
		}
		if( ($exists_opt==1)&&(-e $output_file)){
		   print "\n# write_fasta_seq_by_seq: File $output_file exists, NO write\n";
		}elsif( $out_file_name_given == 1){
		   &write_fasta(\%temp_hash, \$output_file, $rename_seq_opt);
		}else{
		   &write_fasta(\%temp_hash, \$output_file, $rename_seq_opt);
		}
	 }
}
#________________________________________________________________________
# Title     : handle_arguments
# Usage     : Just put the whole box delimited by the two '###..' lines below
#             to inside of your subroutines. It will call 'handle_arguments'
#             subroutine and parse all the given input arguments.
#             To use, claim the arguments, just use the variable in the box.
#             For example, if you had passed 2 file names for files existing
#             in your PWD(or if the string looks like this: xxxx.ext),
#             you can claim them by $file[0], $file[1] in
#             your subroutine.
# Function  : Sorts input arguments going into subroutines and returns default
#             arrays of references for various types (file, dir, hash, array,,,,)
#             If you give (\@out, @file), it will put @out into @array as a ref
#             and also the contents of @out will be dereferenced and put to
#             raw_string regardless what is in it).
#
# Example   : 'handle_arguments(\@array, $string, \%hash, 8, 'any_string')
# Warning   :
# Keywords  : handling arguments, parsing arguments,
# Options   :
# Returns   : Following GLOBAL variables
#
#             $num_opt,    @num_opt     @file          @dir
#             $char_opt,   @char_opt    %vars          @array,
#             @hash        @string,     @raw_string    @range,
#
#             $num_opt has 10,20
#             @num_opt has (10, 20)
#             @file has  xxxx.ext
#             @dir has  dir  or /my/dir
#             $char_opt has 'A,B'
#             @char_opt has (A, B)
#             @array has  (\@ar1, \@ar2)
#             @hash has (\%hash1, \%hash2)
#             @string  ('sdfasf', 'dfsf')
#             @raw_string (file.ext, dir_name, 'strings',,)
#             @range has values like  10-20
#             %vars deals with x=2, y=3 stuff.
#
# Argument  : any type, any amount
# Version   : 4.8
#--------------------------------------------------------------------
sub handle_arguments{
	my($c, $d, $e, $f, $i, $j, $k, $l, $s, $t, $x, $y, $z, $char_opt, $dir, @hash,
		$file, $in_dir, $num_opt, @char_opt, @dir, @file, @string, @file_dir, @k,
		@num_opt, @raw_string,@string, @array, %vars, @range, @temp, $temp,
		@char_options);

  &set_debug_option;
  if(@_<1){ print chr(7),"\n This is handle_arguments. No args Passed, Error?\n"}
  elsif( (@_ ==1)&& (ref($_[0]) eq 'ARRAY') ){ # when there is only 1 argument
	  push(@array, $_[0]);
	  push(@k, $_[0]);
  }elsif( (@_==1)&&( !ref($_[0]) ) ){
	  if(-f $_[0]){ push(@file, $_[0]);   push(@string, $_[0]) }
	  elsif(-d $_[0]){ push(@dir, $_[0]); push(@string, $_[0]) }
	  elsif($_[0]=~/^\d+$/){ push(@num_opt, $_[0]); $num_opt.=$_[0] }
	  elsif($_[0]=~/^\w+$/){ push(@string, $_[0]); }
  }elsif(@_ >=1){ @k = @_ }

  #####______Start of  general argument handling______######
  for($k=0; $k < @k ;$k++){
	  if( !ref($k[$k]) ){
		  if($k[$k]=~ /^[\-]?([a-zA-Z]\d*) {0,5}$/){  push(@char_opt, $1); $char_opt .= "$1\,";
		  }elsif($k[$k]=~ /^\-([a-zA-Z]+)$/){          ## When multiple option is given,
			  @char_options = split(/\,|/, $1);  push(@char_opt, @char_options);
			  $char_opt .= join("\,", @char_options); ## '-' should be used. eg. '-HEGI'
		  }elsif($k[$k]=~ /^(\w+)\=(\S* *)$/){  $vars{$1}=$2;  $vars .= "$1\,";
		  }elsif($k[$k]=~ /^(\-?\d+)$/){ push(@num_opt, $1);  $num_opt .= "$1\,";
		  }elsif($k[$k]=~ /^\d+\.?\d*\-\d+\.?\d*$/){  push(@range,  $k[$k] );
		  }elsif(-f $k[$k]){                          push(@file,   $k[$k] );
		  }elsif(-d $k[$k]){                          push(@dir,    $k[$k] );
		  }elsif($k[$k]=~ /\/[\w\d\.\-]+[\/].+[\/]$/){push(@dir,    $k[$k] );
		  }elsif($k[$k]=~ /^\/[\w\d\.\-]+[\/]*$/){    push(@dir,    $k[$k] );
		  }elsif($k[$k]=~ /^[\/\w\d\-\.]+\.\w+$/){    push(@file,   $k[$k] );
		  }elsif($k[$k]=~ /\S\/[\/\w\d\-\.]+\.\w+$/){ push(@file,   $k[$k] );
		  }elsif($k[$k]=~/^\w+[\/\\\w\d\.\-]+$/){     push(@string, $k[$k] );
		        # string does not have space, but includes '\', '/', '.'
		  }else{                                      push(@raw_string, $k[$k] );  }

	  }elsif( ref($k[$k]) ){
		  if( ref($k[$k]) eq "SCALAR"){
			 if(${$k[$k]} =~ /^[\-]?([a-zA-Z]\d*) {0,5}$/){ push(@char_opt, $1); $char_opt  .= "$1\,";
				}elsif(${$k[$k]}=~ /^\-([a-zA-Z]+)$/){ push(@char_opt, @char_options);
					$char_opt  .= join("\,", @char_options);  ## as an option string.
				}elsif(${$k[$k]}=~ /^(\w+)\=(\S* *)$/){  $vars{$1}=$2;  $vars .= "$1\,";
				}elsif(${$k[$k]}=~ /^(\-?\d+)$/){ $num_opt .= "$1\,";  push(@num_opt, $1);
			    }elsif(${$k[$k]}=~ /^\d+\.?\d*\-\d+\.?\d*$/){    push(@range,  $k[$k] );
				}elsif(-f ${$k[$k]}){                            push(@file,   ${$k[$k]} );
				}elsif(-d ${$k[$k]}){                            push(@dir,    ${$k[$k]} );
				}elsif(${$k[$k]}=~ /\/[\/\w\d\.\-]+[\/].+[\/]/){ push(@dir,    ${$k[$k]} );
				}elsif(${$k[$k]}=~/^\/[\/\w\d\.\-]+[\/]*$/){     push(@dir,    ${$k[$k]} );
				}elsif(${$k[$k]}=~ /^[\/\w\d\-\.]+\.\w+$/){      push(@file,   ${$k[$k]} );
				}elsif(${$k[$k]}=~/^\w+[\w\d\.\-]+$/){           push(@string, ${$k[$k]} );
				}else{                                           push(@raw_string, ${$k[$k]}); }
		  }elsif(ref($k[$k]) eq "ARRAY"){ my @temp_arr = @{$k[$k]}; push(@array, $k[$k]);
			for ($i=0; $i<@temp_arr; $i++){
			   if(-f $temp_arr[$i]){                            push(@file, $temp_arr[$i]);
			   }elsif($temp_arr[$i]=~/^\d+\.?\d*\-\d+\.?\d*$/){ push(@range,$temp_arr[$i] );
			   }elsif(-d $temp_arr[$i]){                        push(@dir , $temp_arr[$i]);
			   }elsif($temp_arr[$i]=~/\/[\/\w\d\.\-]+[\/].+[\/]/){ push(@dir, $temp_arr[$i] );
			   }elsif($temp_arr[$i]=~/^\/[\/\w\d\.\-]+[\/]*$/){ push(@dir, $temp_arr[$i] );
			   }elsif($temp_arr[$i]=~/^[\/\w\d\-\.]+\.\w+$/){   push(@file,$temp_arr[$i] );
																push(@string,$temp_arr[$i] );
			   }elsif($temp_arr[$i]=~/^\w+[\w\d\.\-]+$/){       push(@string,$temp_arr[$i]);
			   }else{                                           push(@raw_string, $temp_arr[$i]); }
			 }
		  }elsif(ref($k[$k]) eq "HASH"){                             push(@hash,   $k[$k] ); }
	  }
  }
  @raw_string=(@raw_string, @string);
  @file = @{&remove_dup_in_arrayH(\@file)};
  #-----------------------------------------------------
	 sub remove_dup_in_arrayH{  my($i, @nondup, @out_ref, %duplicate, @orig, @out_ref);
		for($i=0; $i<@_; $i++){  undef(%duplicate);
	       if(ref($_[$i]) eq 'ARRAY'){    @orig = @{$_[$i]};    }
		   @nondup = grep { ! $duplicate{$_}++ } @orig; push(@out_ref, \@nondup);  }
		if(@out_ref ==1){ return($out_ref[0]);}
		elsif(@out_ref >1){  return(@out_ref);}
	 }
  #-----------------------------------------------------
  return(\@hash, \@array, \@string, \@dir, \@file, \@num_opt,
			\@char_opt, \$num_opt, \$char_opt, \@raw_string, \%vars, \@range );
}
#______________________________________________________________________________
# Title     : sso_to_msp
# Usage     : &sso_to_msp(@ARGV, $single_out_opt);
# Function  : This takes sso file(s) and produces MSP file. It
#             concatenate sso file contents when more than one
#             sso file is given.
# Example   : &sso_to_msp(@ARGV, 'OUT.msp', $single_out_opt);
# Warning   : This capitalize all the input file names when
#              producing xxxxx.msp. xxxxx.sso -> XXXX.sso
# Keywords  : sso_file_to_msp_file, convert_sso_to_msp,
# Options   : _  for debugging.
#             #  for debugging.
#             v  for showing the MSP result to screen
#             s  for making single MSP file for each sso file
#                    as well as big MSP file which has all sso
#             u= for upper expectation value limit
#             l= for lower expect val limit
#             s= for single file name input eg. "s=xxxxx.msp"
#             n  for new format (msp2 format)
#             r  for adding range
#             r2 for adding ranges in all sequence names
#
# Returns   : the file names created (xxxx.msp, yyyy.msp,,,,)
# Argument  :
# Version   : 2.6
#-----------------------------------------------------------------------------
sub sso_to_msp{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
   my ($upper_expect_limit, $lower_expect_limit)=(50, 0);
   my (%sso, @sso, @SSO, $big_out_msp1,  @final_out, $big_out_msp2,
	   $create_sso, $single_out_opt, $add_range, $add_range2, $big_out_msp,
	   $Evalue_thresh, $new_format, $Score_thresh, $margin, $single_file_name);
	if($vars{'u'}=~/([\.\d]+)/){ $upper_expect_limit = $vars{'u'} };
	if($vars{'l'}=~/([\.\d]+)/){ $lower_expect_limit = $vars{'l'} };
	if($vars{'t'}=~/(\d+)/){ $Score_thresh  = $vars{'t'} };
	if($vars{'m'}=~/(\d+)/){ $margin  = $vars{'m'} };
	if($vars{'s'}=~/\S/){ $single_file_name  = $vars{'s'} };
	if($char_opt=~/r2/){  $add_range='r'; $add_range2='r2' }
	if($char_opt=~/r/){   $add_range = 'r' }
	if($char_opt=~/c/){   $create_sso = 'c' }
	if($char_opt=~/s/){   $single_out_opt='s' }
	if($char_opt=~/n/){   $new_format='n' }
   print "\n# File given to sso_to_msp is \"@file\", Normally xxx.sso file names\n";

   if($single_file_name=~/\S/){
	   $big_out_msp=$single_file_name;
   }else{
	   for($i=0; $i < @file; $i++){
		   if($file[$i]=~/\.msp$/){ ## when output file name is given
			   $big_out_msp=$file[$i];
			   splice(@file, $i, 1);
			   $i--;
		   }elsif($file[$i]=~/^(\d+\-\d+)([_\d]*)\.[mfs]?sso/){  ## creates xxxx.msp file name from xxxx.sso
			   $big_out_msp1="\U$1"."$2"."\.msp";
			   $big_out_msp2="\U$1".".msp";
		   }elsif($file[$i]=~/^(\S+)\.[mfs]?sso$/){
			   $big_out_msp1="\U$1"."\.msp";
			   $big_out_msp2="\U$1"."_all".".msp";
			   print "\n# sso_to_msp: File matched  xxxx.sso  format \n";
		   }elsif($file[$i]=~/^(\S+)\.out$/){
			   $big_out_msp1="\U$1"."\.msp";
			   $big_out_msp2="\U$1"."_all".".msp";
			   print "\n# sso_to_msp: File matched  xxxx.out  format \n";
		   }elsif($file[$i]=~/^(\S+)\.p[rot\,]*\.ts\.gz/){
			   $big_out_msp1="\U$1".".msp";
			   $big_out_msp2="\U$1"."_all".".msp";
		   }elsif($file[$i]=~/^(\S+)\.ts\.gz/){
			   $big_out_msp1="\U$1".".msp";
			   $big_out_msp2="\U$1"."_all".".msp";
		   }elsif($file[$i]=~/^(\S+)\.out\.gz/ or $file[$i]=~/^(\S+)\.[mfs]?sso\.gz/){
			   $big_out_msp1="\U$1".".msp";
			   $big_out_msp2="\U$1"."_all".".msp";
		   }
	   }
   }
   if(defined($big_out_msp)){
	   $big_out_msp1=$big_out_msp2=$big_out_msp;
	   print "\n# \$big_out_msp is defined as \'$big_out_msp\'\n";
   }else{
	   print "\n# sso_to_msp: You did not define the big MSP file out format, so $big_out_msp1 \n";
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  (1) When File was given to this sub routine
   #__________________________________________
   if(@file == 1){   ## ONE single file input??
	  print "# one file @file is given, OUT will be: $big_out_msp1 \n";
	  @sso=@{&open_sso_files(@file, $add_range, $add_range2,
	          "u=$upper_expect_limit",
			  "l=$lower_expect_limit",
			  "m=$margin",
			  $new_format,
			  "s=$big_out_msp")};
	  push(@final_out, &write_msp_files(@sso, $big_out_msp1,
	        $single_out_opt, $add_range) );

   }elsif(@file > 1){ ## MOre than 1 file input??
	  @sso=@{&open_sso_files(@file, $add_range, $add_range2,
	        "l=$lower_expect_limit",
	        "u=$upper_expect_limit",
	        "m=$margin",
	        $new_format)};
	  push(@final_out, @{&write_msp_files(@sso, $big_out_msp2,
			$single_out_opt, $add_range)} ); ## concatenates all the hash ref to one
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  (2) When NO File but ARRAY is given
   #      Here, you can have SSO files created
   #__________________________________________
   elsif(@array >=1){
	  print "\n# In sso_to_msp, \@array is given rather than \@file";
	  @sso=@{&open_sso_files(@array, "u=$upper_expect_limit", $add_range2,
			  "l=$lower_expect_limit", $add_range, $create_sso,
			  "m=$margin", $new_format)};
	  push(@final_out, @{&write_msp_files(@sso, $big_out_msp,
						  $single_out_opt, $add_range)} );
   }
   return(\@final_out);
}
#________________________________________________________________________
# Title     : show_options
# Usage     : &show_options;  usually with 'parse_arguments' sub.
# Function  :
# Example   :
# Keywords  : display_options, show_help_options, show_argument_options,
#             show_options_in_headbox, show_prompt_options
# Options   :
# Version   : 1.2
#--------------------------------------------------------------------
sub show_options{
  my($i, @keys, $perl_dir, $arg_num_limit, $head ,$arg_num_limit,
	 @entries_I_want_write );
  my($logname)=getlogin();
  my($pwd)=`pwd`;
  my($date)=`date`;
  chomp($date,$pwd);
  my($not_provided)="--- not provided ---\n";
  my($file_to_read) = $0;

  for($i=0; $i < @_; $i ++){
	  if((ref($_[$i]) eq 'SCALAR')&&(${$_[$i]} =~ /^\d$/)){
		  $arg_num_limit = ${$_[$i]};  }
	  elsif( (!(ref($_[$i]))) && ($_[$i] =~ /^\d$/)){
		  $arg_num_limit = $_[$i];     }
  }
  my %entries = %{&read_head_box(\$file_to_read )};
  if($option_tb_found ==1){
	 @option_tb=@{&read_option_table(\$file_to_read)};
  }
  foreach $help_item (keys %entries){
	  ${$help_item}= $not_provided if( (${$help_item}=~/^[\W]*$/)||( !defined(${$help_item})) );
  }
  #""""""""""""""""""""""""""""""""""""""""
  #########  Writing the format <<<<<<<<<<<
  #""""""""""""""""""""""""""""""""""""""""
  $~ =HEADER_HELP;
  write;   ## <<--  $~ is the selection operator
  $~ =DEFAULT_HELP_FORM;

  @entries_I_want_write=('Options');

  for( @entries_I_want_write ){  write  }

  print chr(7);  print "_"x72,"\n";

  if(@ARGV < $arg_num_limit){ print "\* $0 fataly needs $arg_num_limit arguments\n\n" }

  if(  $option_tb_found == 1){
	 #########  Printing the OPTION table contents <<<<<<<<<<<<
	 print "  Press \"Return\" key to see what options $logname\'s \n\n    \"$0\" take... \n";
		 $key_press=getc();
	 print @option_tb, "\n"x2 if(@option_tb > 0);
  }
format HEADER_HELP  =

**---------------------------------------------------------------------
	O P T I O N S  (I am &show_options)
**--------------------------------------------------------------------
.
format DEFAULT_HELP_FORM =
 @<<<<<<<<<: @*
 $_        $entries{$_}
.
}
#______________________________________________________________
# Title     : write_msp_files
# Usage     : &write_msp_files(\%in1, \%in2, ['s'], [$filename],,)
# Function  : Writes input which is already in msp file format to
#              files either the name is given or generated
#              If more than one ref of hash is given, this will
#              concatenate all the hashes to one big one to
#              make one file.
#             When NO output xxx.msp file name is given, it creates
#              with the query sequence name.
# Example   :  &write_msp_files(@sso, 's', $out_file);
# Warning   : When NO output xxx.msp file name is given, it creates
#              with the query sequence name.
# Keywords  : write_msp,
# Options   : _  for debugging.
#             #  for debugging.
#             s  for each single file output for each hash input
#      filename  for putting output to the specified filename, should be xxx.msp
#
# Returns   : if 's' option is set, it will make say,
#               HI001.msp HI002.msp HI003.msp  rather than
#
#               HI001HI002HI003.msp
#  eg of one output(single file case)
#
#   1027     0.0     1     154   HI0004     1     154   HI0004
#   40       0.0     84    132   HI0004     63    108   HI0001
#   31       0.0     79    84    HI0004     98    103   HI0003
#
# Version   : 2.3
#--------------------------------------------------------------
sub write_msp_files{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

	my ($out_msp_file, $add_range, @final_out, $msp_file_out,
	     @keys, $N, $temp_1, %hash, $query_seq_name, $single_out_opt);

	if($char_opt=~/r/){ $add_range      ='r' };
	if($char_opt=~/s/){ $single_out_opt ='s' };
	if(@file == 1){ $out_msp_file=$file[0]; $single_out_opt='' } # s is for single file output

	if($single_out_opt eq 's'){ #~~~~~~~~~~~` single files output option WITHOUT a given outfilename
		$msp_file_out='default_single_out.msp';
		for($i=0; $i< @hash; $i++){
			my %hash=%{$hash[$i]};
			my @keys =sort keys %hash;
			for($j=0; $j< @keys; $j++){
				#------------------ Writing the first line ---------------------------
				if($keys[$j]=~/(\S+)_\d+\-\d+/){ $N = $1 }else{ $N = $keys[$j] }
				if($hash{$keys[$j]}=~/ +$N[\_\d+\-\d+]* +\d+ +\d+ +$N[\_\d+\-\d+]*/){
				   open(MSP, ">$msp_file_out") ||
					   die "# write_msp_files: I can not create $msp_file_out, check permission\n";
				   chomp( $hash{$keys[$j]} );
				   print MSP $hash{$keys[$j]}, "\n";
				   splice(@keys, $j, 1);
				   $j--; last;
				}
			}
			for($j=0; $j< @keys; $j++){
				chomp( $hash{$keys[$j]} );
				print MSP $hash{$keys[$j]}, "\n";
			}
			print MSP "\n";
		}
		if(-s $msp_file_out){
			print "\n# write_msp_files: $msp_file_out is written \n";
		}else{
			print "\n# Error, write_msp_files\n"; exit
		}
		push(@final_out, $msp_file_out);
		close(MSP);
		return(\@final_out);
	}else{
	   #~~~~~~~~~~~~~ DEfault ~~~~~~~~~~~~~~~~~~
	   #  When output file name was given!
	   #________________________________________
	   if(@file==1){
		   my($temp_1);
		   open(MSP, ">$out_msp_file") ||
			  die "# write_msp_files: I can not create $out_msp_file, check permission\n";
	       for($i=0; $i< @hash; $i++){
			  my %hash=%{$hash[$i]};
			  my @keys =sort keys %hash;
			  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			  # for Final output
			  #_____________________________
			  push(@final_out, $out_msp_file);
			  for($j=0; $j< @keys; $j++){
				  #~~~~~~~ Writing the first line only ~~~~~~~~~~~~~~~~~~
				  if($keys[$j]=~/(\S+)_\d+\-\d+$/){ $N = $1 }else{ $N = $keys[$j] }
				  if($hash{$keys[$j]}=~/ +$N[\_\d+\-\d+]* +\d+ +\d+ +$N[\_\d+\-\d+]*/){
					  $temp_1=$keys[0]; $keys[0]=$keys[$j]; $keys[$j]=$temp_1;
				  }
			  }
			  for($j=0; $j< @keys; $j++){
				  chomp($hash{$keys[$j]});
				  print MSP $hash{$keys[$j]}, "\n";
			  }
			  print MSP "\n";
		   }
		   print MSP "\n";
		   close(MSP);
		   if(-s $out_msp_file and $out_msp_file !~/^ *\.msp$/){
			   print "\n# write_msp_files: $out_msp_file is written\n" if(-s $out_msp_file);
		   }else{
			   print "\n# write_msp_files: ERROR. Either $out_msp_file is empty or \".msp\" is written\n";
		   }
	   }else{
	      for($i=0; $i< @hash; $i++){
			  my %hash=%{$hash[$i]};
			  my @keys =sort keys %hash;
			  ($query_seq_name)=$hash{$keys[0]}=~/\S+ +\d+ +\d+ +(\S+) +\d+ +\d+ +\S+/;
			  $msp_file_out="$query_seq_name\.msp";
			  open(MSP, ">$msp_file_out") or die "\n# write_msp_files: Failed to open $msp_file_out\n";

			  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			  # for Final output
			  #_____________________________
			  push(@final_out, $msp_file_out);
			  for($j=0; $j< @keys; $j++){
				 #~~~~~~~ Writing the first line only ~~~~~~~~~~~~~~~~~~
				 if($keys[$j]=~/(\S+)_\d+\-\d+$/){ $N = $1 }else{ $N = $keys[$j] }
				 if($hash{$keys[$j]}=~/ +$N[\_\d+\-\d+]* +\d+ +\d+ +$N[\_\d+\-\d+]*/){
					$keys[0]=$temp_1; $keys[0]=$keys[$j]; $keys[$j]=$temp_1;
				 }
			  }
			  for($j=0; $j< @keys; $j++){
			     chomp($hash{$keys[$j]});
				 print MSP $hash{$keys[$j]}, "\n";
			  }
			  print MSP "\n";
		   }
		   print MSP "\n";
		   if(-s $out_msp_file and $out_msp_file !~/^ *\.msp$/){
			   print "\n# write_msp_files: $out_msp_file is written\n" if(-s $out_msp_file);
		   }else{
			   print "\n# write_msp_files: ERROR. Either $out_msp_file is empty or only \".msp\" is written\n";
		   }
		   close MSP;
	   }
   }
   if(@final_out > 1){
	   return(\@final_out);
   }else{
	   return(\$final_out[0]);
   }
}
#________________________________________________________________________
# Title     : get_base_names
# Usage     : $base =${&get_base_names(\$file_name)};
#             :   or @bases = &get_base_names(\@files);  # <-- uses `pwd` for abs directory
# Function  : produces the file base name(eg, "evalign"  out of "evalign.pl" ).
# Example   : $base => 'test'  with 'test.txt' or '/home/dir/of/mine/text.txt'
# Warning   :
# Keywords  : get_base_name{, base_name, file_base_name ,  get_file_base_name
#             get_basename, basename, get_root_name
# Options   :
# Returns   :
# Argument  : handles both ref and non-ref.
# Version   : 1.3
#--------------------------------------------------------------------
sub get_base_names{
	my($x, $pos, $pos1, @out_file, $file_only, $file, @file, $base, @base);
	@file=@{$_[0]} if (ref($_[0]) eq 'ARRAY');
	@file=@_ if !(ref($_[0]) eq 'ARRAY');
	for($x=0; $x < @file; $x ++){
		if( ref($file[$x]) ){
			$file = ${$file[$x]};
			$pos1=rindex($file, "/");
	        $file_only=substr($file, ($pos1+1));
			$pos = rindex($file_only, ".");
	        $base= substr($file_only, 0, $pos);
		}else{
			$file = $file[$x];
			$pos1=rindex($file, "/");
	        $file_only=substr($file, ($pos1+1));
			$pos = rindex($file_only, ".");
	        $base= substr($file_only, 0, $pos);
		}
		push(@base, $base);
	}
	if(@base == 1 ){ \$base[0] }else{ \@base }
}
#________________________________________________________________________________
# Title     : read_machine_readable_sso_lines
# Usage     : @out_refs=@{&read_machine_readable_sso_lines(\@SSO, $get_alignment,
#                           $create_sso, $upper_expect_limit,$new_format, $lower_expect_limit,
#                           $attach_range_in_names, $attach_range_in_names2)};
# Function  :
# Example   :
# Keywords  : read_m10_sso_lines, read_msso_lines
# Options   : a c r r2 n
# Version   : 1.2
#--------------------------------------------------------------------------------
sub read_machine_readable_sso_lines{
   my ($upper_expect_limit, $lower_expect_limit)=(50,0);
   my (%match, @out_refs, $target_found, $target_sq_stop, $target_sq_statrt, $match_found,
      $match_seq, $match_found2, $i, $j,$match_found3, $overlap, $sw_score,
      $match_sq_stop, $match_seq2, $sw_ident, $name_range, $target_seq,
      $al_display_start, $match_seq_count);
   for($i=0; $i< @_; $i++){
       if($_[$i]=~/u=(\S+)/){    $upper_expect_limit=$1 }
       elsif(ref($_[$i]) eq 'ARRAY'){ @SSO=@{$_[$i]};   }
       elsif($_[$i]=~/l=(\S+)/){ $lower_expect_limit=$1 }
       elsif($_[$i]=~/^c *$/){     $create_sso = 'c'; print "\n# read_machine_readable_sso_lines: \$create_sso is set"; }
       elsif($_[$i]=~/^a *$/){     $get_alignment='a'; }
       elsif($_[$i]=~/^r *$/){     $attach_range_in_names='r' }
       elsif($_[$i]=~/^r2 *$/){    $attach_range_in_names2='r2' }
       elsif($_[$i]=~/^n *$/){     $new_format='n' }
   }

   print "\n# read_machine_readable_sso_lines : You put PARSEABLE form of sso file";
   for($j=0; $j< @SSO; $j++){
	  if($SSO[$j]=~/\>\>\> *(\S+)\,? +(\d+) +/){  ## >>>  line
		     $target_found=1;  $target_seq_leng=$2;  ## Ignoring the $1, as file name can be different from rea seq names
			 $j+=8;
	  }elsif( $target_found==1 and $SSO[$j]=~/\>\>(\w[\w\-\.]+)([\.prot\,\:]*) */ ){ ##
			 $match_found=1;
			 $match_seq_count++;
			 $al_display_start=0;
			 if(length($2)>0){  print "\n# read_machine_readable_sso_lines: Seq name has this special char \"$2\". I ignore it"; }
			 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			 #  Changing the CASE according to the option
			 #_____________________________________________
			 if($uppercase_seq_name eq 'U'){
				 $match_seq="$1"; $match_seq="\U$match_seq";  ## make it uppercase
			 }elsif($lowercase_seq_name eq 'L'){
				 $match_seq="$1"; $match_seq="\L$match_seq"; ## make it lowercase
			 }else{ $match_seq="$1"; } ## make it uppercase
			 $j+=2;		  next;
	  }elsif($match_found and $SSO[$j]=~/^\; +\w+_expect\:? +(\S+)/){
			 #~~~~~~~~~~~~~~~~~~~~~~~
			 # Filtering by E val
			 #_______________________
			 if( $1 > $upper_expect_limit or $1 < $lower_expect_limit ){
				 $match_found=0; next;
			 }else{ $expect =$1; }
	  }elsif($match_found and $SSO[$j]=~/^ *\; +sw_score *\: +(\S+)/i){  $sw_score =$1;
	  }elsif($match_found and $SSO[$j]=~/^\; +sw_ident\: +(\S+)/){  $sw_ident =$1;
	  }elsif($match_found and $SSO[$j]=~/^ *\; +sw_overlap\: +(\S+)/){  $overlap=$1;
	  }elsif($match_found and $SSO[$j]=~/^ *\>(\w[\w\-\.]+)([\.prot\,\:]*) *[\d+]*/){
			 $match_found2=1;	 $match_found=0;
			 if( length($2)>0 ){  print "\n# read_machine_readable_sso_lines: Seq name has this special char \"$2\". I ignore it"; }
			 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			 #  Changing the CASE according to the option
			 #_____________________________________________
			 if($uppercase_seq_name eq 'U'){
				 $match_seq2="$1"; $match_seq2="\U$match_seq2"; ## make it uppercase
			 }elsif($lowercase_seq_name eq 'L'){
				 $match_seq2="$1";  $match_seq2="\L$match_seq2"; ## make it lowercase
			 }else{ $match_seq2="$1";  }
			 $target_seq=$match_seq2;
	  }elsif($match_found2==1 and $SSO[$j]=~/\; +sq_len\: +(\S+)/){
		     $target_sq_len=$1;
	  }elsif($match_found2==1 and $SSO[$j]=~/\; +al_start\: +(\S+)/){
		     $target_sq_start=$1;
	  }elsif($match_found2==1 and $SSO[$j]=~/\; +al_stop\: +(\S+)/){
		     $target_sq_stop=$1;
	  }elsif($SSO[$j]=~/\; +al_display_start/ and $al_display_start < 1){
             $al_display_start ++;
	  #------------------------------------------------------------
	  }elsif($match_found2 and $SSO[$j]=~/\>(\w[\w\-\.]+)([\.prot\,\:]*) *[\d+]*/){
             $match_found3=1; $match_found2=0;
             if(length($2)>0){  print "\n# open_sso_files: Seq name has this special char \"$2\". I ignore it"; }
	  }elsif($match_found3 and $SSO[$j]=~/\; +sq_len\: +(\S+)/){
		     $match_sq_len=$1;
	  }elsif($match_found3 and $SSO[$j]=~/\; +al_start\: +(\d+)/){
		     $match_sq_start=$1;
	  }elsif($match_found3 and $SSO[$j]=~/\; +al_stop\: +(\d+)/){
             $match_sq_stop=$1;
	  }elsif($match_found3 and $SSO[$j]=~/\; +al_display_start/){
			 $match_found3=0;          $al_display_start++;
			 if($expect=~/^$/){ $expect='0.0'; }
			 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			 # adding the offset for names with ranges
			 #__________________________________________________
			 if($target_seq=~/^\S+_(\d+)\-(\d+)/){ $target_sq_start +=$1-1; $target_sq_stop +=$1-1;  }

             #~~~~~~~~~~~~~~~~~~~~~~~~~
             # Attaching the ranges  (under NO e option)
             #_________________________
             if($attach_range_in_names==1){
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  # Checks margin opt and adds it
                  #__________________________________
                  if($margin=~/\d+/){
                      if($match_sq_start < $margin){  $match_sq_start=1;
                      }else{          $match_sq_start-=$margin;   }
                      $match_sq_stop += $margin;
                  }
                  $name_range="$match_seq\_$match_sq_start\-$match_sq_stop";

                  #~~~~~~~~ If 'rr' opt is set, put ranges for both target and match seqs ~~~~~~~
                  if($attach_range_in_names2==1 and $target_seq !~/^\S+_(\d+)\-(\d+)/){
                      $target_seq="$target_seq\_$target_sq_start\-$target_sq_stop";
                  }
                  if($original_target_seq=~/\S+/){ $target_seq=$original_target_seq } # for PVM version out
                  if($new_format=~/n/){  # under NO e option
                      $match{$name_range}=
                         sprintf("%s %s %s %s %s %s %s %s %s\n",
                         $target_seq, $target_sq_start, $target_sq_stop, $sw_score, $expect, $sw_ident,
                         $match_sq_start, $match_sq_stop, $name_range);
                  }else{
                      $match{$name_range}=
                         sprintf("%-5s %-8s %-6s %-4s %-5s %-30s %-4s %-5s %s\n",
                         $sw_score, $expect, $sw_ident, $target_sq_start, $target_sq_stop, $target_seq,
                         $match_sq_start, $match_sq_stop, $name_range);
                  }
             }else{
                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 # Checks margin opt and adds it
                 #__________________________________
                 if($margin=~/\d+/){
                     if($match_sq_start < $margin){  $match_sq_start=1;
                     }else{                          $match_sq_start-=$margin; }
                     $match_sq_stop += $margin;
                 }
                 if($original_target_seq=~/\S+/){ $target_seq=$original_target_seq } # for PVM version out
                 if($new_format=~/n/){
                     $match{$match_seq}=
                        sprintf("%s %s %s %s %s %s %s %s %s\n",
                        $target_seq, $target_sq_start, $target_sq_stop, $sw_score, $expect, $sw_ident,
                        $match_sq_start, $match_sq_stop, $match_seq);
                 }else{
                    $match{$match_seq}=sprintf("%-5s %-8s %-6s %-4s %-5s %-30s %-4s %-5s %s\n",
                       $sw_score, $expect, $sw_ident, $target_sq_start, $target_sq_stop, $target_seq,
                       $match_sq_start, $match_sq_stop, $match_seq);
                 }
             }
	  }elsif($get_alignment and $al_display_start==1 and $SSO[$j]=~/^([\w\-]+) *$/){
		  ${"match_alignment\_$match_seq_count"}{$match_seq2} .= $1;
	  }elsif($get_alignment and $al_display_start==2 and $SSO[$j]=~/^([\w\-]+) *$/){
		  ${"match_alignment\_$match_seq_count"}{"$match_seq"} .= $1;
	  }elsif($get_alignment and $SSO[$j]=~/^ *\;al_cons\:/){
		  $al_display_start=0;
		  my %temp=%{"match_alignment\_$match_seq_count"};
		  push(@out_refs, \%temp );
		  %{"match_alignment\_$match_seq_count"}=();
	  }
   } ## <-- for($j=0; $j< @SSO; $j++)

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # If create sso option is set, it creates SSO files(array input case)
   #________________________________________________________________________
   if( $create_sso and !$get_alignment){
	   open (SSO2, ">$target_seq\.msso");
	   print SSO2 @SSO, "\n";
       print "\n# read_machine_readable_sso_lines : $target_seq\.msso file  is created by \"c\" opt ";
	   close SSO2
   }
   unless($get_alignment){
	   push(@out_refs, \%match);
   }
   return(\@out_refs);
}

#________________________________________________________________________
# Title     : show_hash
# Usage     : &show_hash(\@input_array);
# Function  : for debugging purpose. Shows any array elem line by line.
#             the line is 60 elements long (uses recursion)
# Example   : Output:      item1
#             Output:      item2
#             Output:      item3
# Warning   : There is a global variable:  $show_hash_option
#             It tries to detect any given sting which is joined by ','
# Keywords  :
# Options   : -s or -S or s or S for spaced output. Eg)
#             seq1       1 1 1 1 1 1 1 1 1 1 1 1
#
#             instead of
#             seq1       111111111111
#
#             -h or -H or h or H for horizontal line of '---------...'
#
# Returns   :
# Argument  :
# Version   : 1.7
#--------------------------------------------------------------------
sub show_hash{
  my($k, $i, $t, @in2, $in, $LEN, %TEM ); ## You should not put $show_hash_option
  my(@in)=@_;                     ## and $horizontal_line  in my !!!
  my($KL)=2; # default keys string length;
  my($VL)=80; # default values string length;
  my($GAP)=2;  # default space between keys and values
  my($horizontal_line, $show_hash_optionXX, $Hash_counter, @line);

  ## This is to get the option of 'space' to make spaced output.
  for($t=0; $t < @in; $t++){
	 if($in[$t] =~/^[-]+[sS][pace]*$/){
		 $show_hash_optionXX = 1;
		 splice(@in, $t, 1);
	 }elsif(${in[$t]} =~/^[-]+[sS][pace]*$/){
		 $show_hash_optionXX = 1;
		 splice(@in, $t, 1);
	 }elsif($in[$t] =~/^[-]+[hH][rR]*$/){
		 $horizontal_line = 1;
		 splice(@in, $t, 1);
	 }elsif(${in[$t]} =~/^[-]+[hH][rR]*$/){
		 $horizontal_line = 1;
		 splice(@in, $t, 1);
	 }
  }

  ######## Main loop #################
  if($horizontal_line ==1){  ## This puts the delimiter '--------------(  )'
	  $Hash_counter ++;
	  print "\n","-"x78,"(${Hash_counter}th hash)", "\n";
  }

  for($k=0; $k < @in; $k++){
	 if(ref($in[$k]) eq 'ARRAY'){  ### When the hashes were given in array ref.
		 &show_hash(@{$in[$k]}, $show_hash_optionXX, $horizontal_line);
		 print "\n";
	 }
	 elsif(ref($in[$k]) eq 'HASH'){  ### recursion
		 &show_hash(%{$in[$k]});
		 print "\n";
	 }
	 elsif(ref($in[$k+1]) eq 'HASH'){  ### recursion
		 &show_hash(%{$in[$k+1]}); print "\n";
	 }
	 elsif(ref($in[$k]) eq 'SCALAR'){  print ${$_[$k]}, "\n";  }
	 elsif( !ref($in[$k]) ){
		 if( !ref($in[$k+1]) && defined($in[$k+1])  ){
			 if($show_hash_optionXX == 1){  #### space option checking.
				#if($in[$k+1] =~ /\,.+\,/){  #### if the string is joined with ','
				#	 @line = split(/\,/, $_[$k+1]);
				# }else{
				#	 @line = split(//, $_[$k+1]);
				# }
				%TEM = @in;
				$LEN = ${&max_elem_string_array_show_hash(keys %TEM)};
				 if($LEN > $KL){ $KL = $LEN + $GAP +2};
				 printf ("%-${KL}s ", $in[$k]);  $k++;
				 printf ("%-${VL}s\n","@line");
			 }else{                        ### If not option is set, just write
				%TEM = @in;
				$LEN = ${&max_elem_string_array_show_hash( keys %TEM)};
				 if($LEN > $KL){ $KL = $LEN + $GAP +2};
				 printf ("%-${KL}s ", $in[$k]);  $k++; # print $in[$k], "\t";  $k++;
				 printf ("%-${VL}s\n",$in[$k]);        # print $in[$k], "\n";
			 }
		 }
		  #________________________________________________________
		  # Title    : max_elem_string_array_show_hash
		  # Keywords : largest string length of array
		  # Function : gets the largest string length of element of any array of numbers.
		  # Usage    : ($out1, $out2)=@{&max_elem_array(\@array1, \@array2)};
		  #            ($out1)       =${&max_elem_array(\@array1)          };
		  # Argument : numerical arrays
		  # returns  : one or more ref. for scalar numbers.
		  # Version  : 1.1
		  #-------------------------------------------------------
		  sub max_elem_string_array_show_hash{
			 my(@input, $i, $max_elem);
			 @input = @{$_[0]} || @_;
			 for($i=0; $i< @input ; $i++){
					$max_elem = length($input[0]);
					if (length($input[$i]) > $max_elem){
						$max_elem = length($input[$i]);
					}
			 }
			 \$max_elem;
		  }
		  #####################################insert_gaps_in_seq_hash
	 }
  }
}



#________________________________________________________________________________
# Title     : read_machine_unreadable_sso_lines
# Usage     : @out_refs=@{&read_machine_unreadable_sso_lines(\@SSO, $get_alignment,
#                           $create_sso, $upper_expect_limit,$new_format, $lower_expect_limit,
#                           $attach_range_in_names, $attach_range_in_names2)};
# Function  :
# Example   :
# Keywords  : read_normal_sso_lines
# Options   : a c r r2 n
# Version   : 1.0
#--------------------------------------------------------------------------------
sub read_machine_unreadable_sso_lines{
   my ($upper_expect_limit, $lower_expect_limit)=(50,0);
   my (@SSO, @out_refs, $match_seq, $match_evalue, $target_found, $target_seq_len, $space, %match);
	  for($i=0; $i< @_; $i++){
		  if($_[$i]=~/u=(\S+)/){    $upper_expect_limit=$1 }
		  elsif(ref($_[$i]) eq 'ARRAY'){ @SSO=@{$_[$i]};   }
		  elsif($_[$i]=~/l=(\S+)/){ $lower_expect_limit=$1 }
		  elsif($_[$i]=~/^c$/){     $create_sso = 'c' }
		  elsif($_[$i]=~/^a$/){     $get_alignment='a'; }
		  elsif($_[$i]=~/^r$/){   $attach_range_in_names='r' }
		  elsif($_[$i]=~/^r2$/){   $attach_range_in_names2='r2' }
		  elsif($_[$i]=~/^n$/){   $new_format='n' }
	  }

   print "\n# open_sso_files : You have put non-parseable format of xxxx.sso\n";
   print "\n#                : Did you set \'M\' option in do_sequence_search? \n";

   for($j=4; $j< @SSO; $j++){
	   if($SSO[$j]=~/^ *\S+\: +(\d+) +\w+/){
		  $target_seq_len=$1;
		  print "\n target seq len is  $target_seq_len \n";
				 # matching  >MJ0497
	   }elsif($SSO[$j]=~/^ \>(\w[\w\-\.\/\\]+)/){
		   $target_seq_name=$1;
		   $j+=3; ## jumping to skip the stat bars
		   print "\n# open_sso_files : Found Query seq=> $target_seq_name ";
				  # matching >>MG032 ATP-d (addA) Bacillus subtilis  (666 aa)
	   }elsif($SSO[$j]=~/^ {0,4}\>\> *(\S+) +.+\((\d+) aa\)$/){
		  $entry_found=1;     $target_found=0;
		  $target_gap_len=0;  $match_gap_len=0;
		  $match_seq=$1;      $match{$match_seq} ="$target_seq_name $target_seq_len $match_seq $2 ";
		  print "\n# open_sso_files : Found MATCHed seq $match_seq\n";
	   }elsif($SSO[$j]=~/expect\( *\) +(\S+)/){ ## getting Evalue
		  $match_evalue=$1;   $match{$match_seq} .="$match_evalue ";
	   }elsif($SSO[$j]=~/Smith\-Waterman +score\: +(\d+)\;.+in (\d+) aa overlap/i){
		  $sw_score=$1;       $overlap=$2;
		  $match{$match_seq}.="$sw_score $overlap ";
	   }elsif( $target_found < 1 and $SSO[$j]=~/^( +)(\d+) +\d+/  ){
		  $gap_start=length($1)+length($2);
		  $start=$2;          $target_found=1;
	   }elsif( $target_found==1 and $SSO[$j]=~/^( +)[\.\:]/ ){ ### matching    .: .: : ::     :.:..: :.. .. ..
		  $space=length($1);
		  $target_seg_start=$space-$gap_start+$start;
		  $target_seg_end=$target_seg_start+$overlap;
		  $target_range="$target_seg_start-$target_seg_end";
	   }elsif( defined($space) and $target_found ==1 and  $SSO[$j]=~/^( +)(\d+)/ ){
		  $target_found++;
		  $match_gap_start=length($1)+length($2);
		  $match_start=$2;
		  $match_seg_start=$space-$match_gap_start+$match_start;
		  $match_seg_end=$match_seg_start+$overlap;
		  $match_range ="$match_seg_start-$match_seg_end";
		  $match{$match_seq}.="$target_range $match_range ";
		  #print "\n $target_seq_name $match_seq $match_evalue $overlap $target_range $match_range";
	   }
	}# end of for $j
	if( ($create_sso=~/c/) && (@file < 1) ){
	   open (SSO2, ">$target_seq_name\.sso");
	   print SSO2 @SSO, "\n";
	   print "\n# $target_seq_name\.sso  is created";
	   close SSO2;
	}
	push(@out_refs, \%match);
	return(\@out_refs);
}# end of for $i
