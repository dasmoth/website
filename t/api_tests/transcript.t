#!/usr/bin/env perl

# This is a unit test template for implementing tests that work
# with a running WormBase Website instance.
#
# Unit tests are called automagically, just adhere to the following:
#
# 1. the unit test is placed in the t/api_tests folder
# 2. the filename and package name coincide (sans suffix)
# 3. unit test names have the prefix "test_"
#
# Actual tests are implemented at the bottom of this file. Please see:
#
# 1. test_single_gene

{
    # Package name is the same as the filename (sans suffix, i.e. no .t ending)
    package transcript;

    # Limit the use of unsafe Perl constructs.
    use strict;

    # We use Test::More for all tests, so include that here.
    use Test::More;

    # This variable will hold a reference to a WormBase API object.
    my $api;

    # A setter method for passing on a WormBase API object from t/api.t to
    # the subs of this package.
    sub config {
        $api = $_[0];
    }

    sub test_transcript_features {
        my $transcript = $api->fetch({ class => 'Transcript', name => 'B0336.10.2' });

        can_ok('WormBase::API::Object::Transcript', ('feature'));

        my $features = $transcript->feature();

        # Please keep test names/descriptions all lower case.
        isnt($features->{'data'}, undef, 'data returned');
        is  ( scalar keys %{$features->{'data'}}, 3, "correct number of features returned");
    }

    sub test_print_sequence {

        # -ve strand test
        my $transcript = $api->fetch({ class => 'Transcript', name => 'B0041.11.1' });
        can_ok('WormBase::API::Object::Transcript', ('print_sequence'));

        my $sequence = $transcript->print_sequence();
        isnt($sequence->{'data'}, undef, 'data returned');
        is (scalar @{$sequence->{'data'}}, 3, 'correct amount of sequences returned -ve strand');

        foreach my $seq (@{$sequence->{'data'}}) {
            if($seq->{'header'} eq "unspliced + UTR"){
                is($seq->{'length'}, 1059, 'correct unspliced length returned -ve strand');
            } elsif ($seq->{'header'} eq "spliced + UTR") {
                is($seq->{'length'}, 427, 'correct spliced length returned -ve strand');
                my $spliced = '<span style="Color: #808080">aggatcggatgg</span><span style="Background-color: #FFFF00">ATGGATGGATCAACGACGATCGGGATGGATAGATCCGGGACAGATCGAGATGGACTG</span><span style="Background-color: #FFA500">ATCTCTAGCAA' . "\n" .
                    'ACGACGAAGGATGGATGGTAGGAAAGGATGGATGGTTGGAACAATGGAAATAGAAAATCCAAAGGG</span><span style="Background-color: #FFFF00">ATGTTGA</span><span style="Color: #808080">aagctga' . "\n" .
                    'ttgaagaaattctggacaatcggttgaatgggccgaccaatcaggatggatggatcgggacgagacgacgacgacgacaa' . "\n" .
                    'acaataaggatggacggatggattgggacgatcgaaagggatttgggacaaaatgacgacgacaaaggatggatgaattc' . "\n" .
                    'gggacggatcgggatggacg</span><span style="Color: #808080">cgatcgatggatggatggacgaggcaacaacgaccacgaacgattaggatggatggatgg' . "\n" .
                    'aaggatcgatcg</span><span style="Color: #808080">gatctcaagcaaaca</span>';
                is($seq->{'sequence'}, $spliced, "correct spliced sequence returned -ve strand");
            } elsif ($seq->{'header'} eq "conceptual translation") {
                is($seq->{'length'}, 46, 'correct protein length returned -ve strand');
                is($seq->{'sequence'}, "MDGSTTIGMDRSGTDRDGLISSKRRRMDGRKGWMVGTMEIENPKGC*", 'correct protein sequence returned -ve strand');
            }
        }

        # +ve strand test
        $transcript = $api->fetch({ class => 'Transcript', name => 'AC3.1' });

        $sequence = $transcript->print_sequence();
        isnt($sequence->{'data'}, undef, 'data returned');
        is (scalar @{$sequence->{'data'}}, 3, 'correct amount of sequences returned +ve strand');

        foreach my $seq (@{$sequence->{'data'}}) {
            if($seq->{'header'} eq "unspliced + UTR"){
                is($seq->{'length'}, 1313, 'correct unspliced length returned +ve strand');
            } elsif ($seq->{'header'} eq "spliced + UTR") {
                is($seq->{'length'}, 1068, 'correct spliced length returned +ve strand');
                my $spliced = '<span style="Background-color: #FFFF00">ATGATCATGTTCACAGAAGCTGAAGTTATGAGTTTTTCATACGCCGTTGATTTTGGAGTTCCCGAATGGCTCAAACTTTA' . "\n" .
                    'CTATCACGTCATTTCCGTGGTGTCAACTGTTATTTCATTTTTCTCAATGTACATAATTTTGTTTCAAAGTGGGAAAATGG' . "\n" .
                    'ATGGATATCGATTCTATCTATTTTATATGCAG</span><span style="Background-color: #FFA500">TTTGCTGGATGGTTGATGGATCTTCATCTATCTACTTTTATGCAGTTC' . "\n" .
                    'ATTCCATTATTCCCAGTTTTTGGAGGATATTGTACTGGACTCTTGACTCAAATTTTCAGAATTGACGATTCATTTCAAAC' . "\n" .
                    'G</span><span style="Background-color: #FFFF00">ACATATACTGCATTTACCATTTGCTTGGTAGCCAGTGCTTTGAATAGTTGCTTTGTTCGGAAGCATCAAGCAATTTCTA' . "\n" .
                    'AAATCAGCTCTAAATATCTTCTCGATAATGTTACATACTGCATTGTTATATTTCTACTCAATATATATCCAGTTATTGCT' . "\n" .
                    'GCATCTTTACTTTATTTGAGCATGCTCAATAAGTCCGAACAAGTTGAATTGGTGAAATCG</span><span style="Background-color: #FFA500">GTTTACCCAAATCTCGTTGA' . "\n" .
                    'TAAATTTGCAAGTCTACCAAACTACGTGGTATTTGATTCCAATATATGGGCAATTGTATTCTTTGCATTCATATTTTTTG' . "\n" .
                    'GTTGTACATATACACTTGTTTTGATTGTCACAACTACTTATCAAATGTTCAAAATATTAGACGATAATCGAAAACATATC' . "\n" .
                    'AGTGCTTCAAACTATGCAAAGCATCGAGCCACTTTGAGAAGTCTTTTAGCTCAGTTTACAACGTGTTTCTTGATTGTTGG' . "\n" .
                    'TCCAGCGTCTTTGTTCTCTTTACTGGTAGTTATAAGATATGAACATAGTCAAG</span><span style="Background-color: #FFFF00">TAGCAACACATTGGACCATTGTTGCTC' . "\n" .
                    'TAACTCTCCATTCCAGTGCAAATGCAATTGTAATGGTTATCACATATCCTCCATACAGACATTTTGTAATGCTATGGAAA' . "\n" .
                    'ACGAACAG</span><span style="Background-color: #FFA500">ATCATTCCACTTCGCATCATCTCAATATCAACGGTCTACTCTCCCGAATACAAGAATTCAAACGGAGCGAAG' . "\n" .
                    'TATTGCAGTAACAATAACAACCCATTAA</span>';
                is($seq->{'sequence'}, $spliced, "correct spliced sequence returned +ve strand");
            } elsif ($seq->{'header'} eq "conceptual translation") {
                is($seq->{'length'}, 355, 'correct protein length returned +ve strand');
                is($seq->{'sequence'}, "MIMFTEAEVMSFSYAVDFGVPEWLKLYYHVISVVSTVISFFSMYIILFQSGKMDGYRFYLFYMQFAGWLMDLHLSTFMQFIPLFPVFGGYCTGLLTQIFRIDDSFQTTYTAFTICLVASALNSCFVRKHQAISKISSKYLLDNVTYCIVIFLLNIYPVIAASLLYLSMLNKSEQVELVKSVYPNLVDKFASLPNYVVFDSNIWAIVFFAFIFFGCTYTLVLIVTTTYQMFKILDDNRKHISASNYAKHRATLRSLLAQFTTCFLIVGPASLFSLLVVIRYEHSQVATHWTIVALTLHSSANAIVMVITYPPYRHFVMLWKTNRSFHFASSQYQRSTLPNTRIQTERSIAVTITTH*", 'correct protein sequence returned +ve strand');
            }
        }


        # non-coding test
        $transcript = $api->fetch({ class => 'Transcript', name => 'B0379.1b' });

        $sequence = $transcript->print_sequence();
        isnt($sequence->{'data'}, undef, 'data returned');
        is (scalar @{$sequence->{'data'}}, 1, 'correct amount of sequences returned nc transcript');
        is ($sequence->{'data'}[0]->{'sequence'}, "atgggcgccttcgatcgcgtgaaagctcaagttgcatccgactcaaaatggacatcagctccttacaagggatttgtggccggaagcccatcaaacacgtatattgatattgtttccactgcgtgagttttcaacatgcaacctatcctgaagctttttaaaattaattctttgaagacgccactaacacgatgaattttgctcgtcactcgaaatacgatgaaatgtattctccttatctcggatcattccgcgaacgacacaattatacttcaattgctccaagcttgtgtattaacaaaacaaaccgtgccatcgagtatgacctggcaccacacaaggcttacaatccacgacaatccgaatggcttcttgaaaaagacaagaaatatagagttcgtggtgctcgtaatttaatttacacaaaaagcgcatcggatatcagtttgcctccactgacacgtcgcacattcacagttccaacagatactcttcgtcatcagaatcaatttctctactggaatggtcgtgcacttggtcttgactatgttgctccattccttcgtcgtgaagattattctcgtcacgaggatcgccgttatcagagaatttactggtctccacatttcattgatttgcttccatcttgccgtcattctgcacatcttatgctttccgcttattaa", 'correct sequence returned nc transcript' )
    }

    #related to issue #2710
    sub test_expression_widget {
        my $transcript = $api->fetch({ class => 'Transcript', name => 'B0336.6.1' });

        #test microarray_topology_map_position
        can_ok('WormBase::API::Object::Transcript', ('microarray_topology_map_position'));

        my $profiles = $transcript->microarray_topology_map_position();
        is($profiles->{'description'}, 'microarray topography map', 'correct description returned');
        isnt($profiles->{'data'}, undef, 'data returned');

        my $p = shift @{$profiles->{'data'}};

        is($p->{'id'}, 'B0336.6', 'correct profile id returned');
        is($p->{'class'}, 'expr_profile', 'correct class returned');
        ok($p->{'label'} =~ /Mountain: 11/, 'correct mountain for expression profile returned');


        #test anatomic_expression_patterns
        can_ok('WormBase::API::Object::Transcript', ('anatomic_expression_patterns'));

        my $patterns = $transcript->anatomic_expression_patterns();
        isnt($patterns->{'data'}, undef, 'data returned');
        is($patterns->{'description'}, 'expression patterns for the gene' , 'correct description returned ');
        is($patterns->{'data'}->{'image'}, '/img-static/virtualworm/Gene_Expr_Renders/WBGene00015146.jpg' , 'correct image returned');


        #test expression_patterns
        can_ok('WormBase::API::Object::Transcript', ('expression_patterns'));
        my $expressions = $transcript->expression_patterns();
        isnt($expressions->{'data'}, undef, 'data returned');
        is($expressions->{'description'}, 'expression patterns associated with the gene:WBGene00015146' , 'correct description returned ');
        is($expressions->{'data'}[0]->{'description'}->{'text'}, 'Collectively, these approaches revealed that ABI-1 is expressed in a number of neurons within the nerve ring and head, including the amphid interneurons AIYL/R, the RMEL/R motoneurons, coelomocytes, and several classes of ventral cord motoneuron.' , 'correct expression description returned' );
        is($expressions->{'data'}[0]->{'type'}, 'Reporter gene', 'correct type returned');
        is($expressions->{'data'}[0]->{'expression_pattern'}->{'id'}, 'Expr8549', 'correct expression pattern returned');

        #test expression_profiling_graphs
        can_ok('WormBase::API::Object::Transcript', ('expression_profiling_graphs'));
        my $graphs = $transcript->expression_profiling_graphs();
        isnt($graphs->{'data'}, undef, 'data returned');
        is($graphs->{'description'}, 'expression patterns associated with the gene:WBGene00015146' , 'correct description returned ');
        is($graphs->{'data'}[0]->{'description'}->{'text'}, 'Developmental gene expression time-course.  Raw data can be downloaded from ftp://caltech.wormbase.org/pub/wormbase/datasets-published/levin2012' , 'correct expression description returned' );
        is($graphs->{'data'}[0]->{'type'}, 'Microarray', 'correct type returned');
        is($graphs->{'data'}[0]->{'expression_pattern'}->{'id'}, 'Expr1011958', 'correct expression pattern returned');

        #test anatomy_terms
        can_ok('WormBase::API::Object::Transcript', ('anatomy_terms'));
        my $anatomy_terms = $transcript->anatomy_terms();
        isnt($anatomy_terms->{'data'}, undef, 'data returned');
        is($anatomy_terms->{'description'}, 'anatomy terms from expression patterns for the gene' , 'correct description returned ');
        is($anatomy_terms->{'data'}->{'WBbt:0005751'}->{'class'}, 'anatomy_term' , 'correct anatomy term class returned');
        is($anatomy_terms->{'data'}->{'WBbt:0005751'}->{'label'}, 'coelomocyte' , 'correct anatomy term label returned');
        is($anatomy_terms->{'data'}->{'WBbt:0005751'}->{'id'}, 'WBbt:0005751' , 'correct anatomy term id returned');

        #test expression_cluster
        can_ok('WormBase::API::Object::Transcript', ('expression_cluster'));
        my $expression_cluster = $transcript->expression_cluster();
        isnt($expression_cluster->{'data'}, undef, 'data returned');
        is($expression_cluster->{'description'}, 'expression cluster data' , 'correct description returned ');
        is($expression_cluster->{'data'}[0]->{'expression_cluster'}->{'id'}, 'cgc4489_group_2' , 'correct expression cluster id returned');
        is($expression_cluster->{'data'}[0]->{'description'}, 'Genome-wide analysis of developmental and sex-regulated gene expression profile.' , 'correct expression cluster description returned');

        #test fpkm_expression_summary_ls
        can_ok('WormBase::API::Object::Transcript', ('fpkm_expression_summary_ls'));
        my $fpkm_expression_summary_ls = $transcript->fpkm_expression_summary_ls();
        isnt($fpkm_expression_summary_ls->{'data'}, undef, 'data returned');
        is($fpkm_expression_summary_ls->{'description'}, 'Fragments Per Kilobase of transcript per Million mapped reads (FPKM) expression data' , 'correct description returned ');
        is($fpkm_expression_summary_ls->{'data'}->{'plot'}, '/img-static/rplots/WS244/4876/fpkm_WBGene00015146.png' , 'correct plot returned');
        is($fpkm_expression_summary_ls->{'data'}->{'table'}->{'fpkm'}->{'data'}[0]->{'project'}, 'Thomas Male Female comparison' , 'correct project returned');



    }

}

1;
