// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ClassicAbinitio.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange
/// @author James Thompson
/// @author Mike Tyka

// Unit Headers
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/simple_moves/SymmetricFragmentMover.hh>

// Package Headers
#include <protocols/simple_moves/GunnCost.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/WhileMover.hh>
#include <protocols/abinitio/AllResiduesChanged.hh>

//for cenrot
#include <protocols/moves/CompositionMover.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
//#include <protocols/simple_moves/BackboneMover.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/ozstream.hh>
#include <numeric/numeric.functions.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>

//000000000000000000000000000000000000000000000000000000000
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <stack>
#include <numeric>
#include <core/scoring/rms_util.hh>
#include <core/import_pose/import_pose.hh> 
#include <core/util/SwitchResidueTypeSet.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
//======================================================
#include <stdio.h>
#include <string.h>
#include <vector>
#include <iterator>
#include <math.h>
// #include <tuple>
#include <string>
//======================================================

using namespace std;
core::pose::Pose initPose;
core::pose::Pose tempPose;
core::pose::Pose native_pose;
//000000000000000000000000000000000000000000000000000000000

//000000000000000000000000000000000000000000000000000000000
///@brief 控制输出格式
#include<iomanip>
//000000000000000000000000000000000000000000000000000000000
//// C++ headers
#include <cstdlib>
#include <string>
#ifdef WIN32
#include <ctime>
#endif

//debug

#include <protocols/moves/MonteCarlo.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


static basic::Tracer tr( "protocols.abinitio" );

using core::Real;
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

/*!
@detail call this:
ClassicAbinitio::register_options() before devel::init().
Derived classes that overload this function should also call Parent::register_options()
*/

// This method of adding options with macros is a pain in the ass for people
// trying to nest ClassicAbinitio as part of other protocols. If you don't call
// ClassicAbinitio::register_options() in your main function, you get a really
// unintuitive segfault as the options system doesn't know about the options
// listed below. The solution is to call register_options() in your main method
// before devel::init(), which is really ugly as the main method shouldn't need
// to know what protocols are called, and it's prone to error because it's an
// easy thing to forget.
// This should get some more thought before it becomes the standard way to add options.

void protocols::abinitio::ClassicAbinitio::register_options() {
	Parent::register_options();
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant( OptionKeys::abinitio::increase_cycles );
	option.add_relevant( OptionKeys::abinitio::smooth_cycles_only );
	option.add_relevant( OptionKeys::abinitio::debug );
	option.add_relevant( OptionKeys::abinitio::skip_convergence_check );
	option.add_relevant( OptionKeys::abinitio::log_frags );
	option.add_relevant( OptionKeys::abinitio::only_stage1 );
	option.add_relevant( OptionKeys::abinitio::end_bias );
	option.add_relevant( OptionKeys::abinitio::symmetry_residue );
	option.add_relevant( OptionKeys::abinitio::vdw_weight_stage1 );
	option.add_relevant( OptionKeys::abinitio::override_vdw_all_stages );
	option.add_relevant( OptionKeys::abinitio::recover_low_in_stages );
	option.add_relevant( OptionKeys::abinitio::close_chbrk );
}


namespace protocols {
namespace abinitio {
	double ClassicAbinitio::TrialEnergy_ = 0.0;

//little helper function
bool contains_stageid( utility::vector1< ClassicAbinitio::StageID > vec, ClassicAbinitio::StageID query ) {
	return find( vec.begin(), vec.end(), query) != vec.end();
}

/// @detail  large (stage1/stage2)
/// small(stage2/stage3/stage4)
/// smooth_small ( stage3/stage4)
ClassicAbinitio::ClassicAbinitio(
	simple_moves::FragmentMoverOP brute_move_small,
	simple_moves::FragmentMoverOP brute_move_large,
	simple_moves::FragmentMoverOP smooth_move_small,
	int  /*dummy otherwise the two constructors are ambiguous */
) :
	brute_move_small_(std::move( brute_move_small )),
	brute_move_large_( brute_move_large ),
	smooth_move_small_(std::move( smooth_move_small ))
{
	BaseClass::type( "ClassicAbinitio" );
	get_checkpoints().set_type("ClassicAbinitio");
	// std::cerr << "ClassicAbinitio::constructor has stubbed out...(fatal) see code file";
	// runtime_assert( 0 ); //---> needs more implementation to use this constructor: e.g. read out movemap from FragmentMover...
	movemap_ = brute_move_large->movemap();
	//  set_defaults( pose ); in constructor virtual functions are not called
	bSkipStage1_ = false;
	bSkipStage2_ = false;

	close_chbrk_ = false;

	stage4_cycles_pack_rate_ = 0.25;
}

ClassicAbinitio::ClassicAbinitio(
	core::fragment::FragSetCOP fragset_small,
	core::fragment::FragSetCOP fragset_large,
	core::kinematics::MoveMapCOP movemap
)  :
	movemap_( movemap )
{
	BaseClass::type( "ClassicAbinitio" );
	get_checkpoints().set_type("ClassicAbinitio");
	using namespace basic::options;
	using simple_moves::FragmentCostOP;
	using simple_moves::ClassicFragmentMover;
	using simple_moves::SymmetricFragmentMover;
	using simple_moves::SmoothFragmentMover;
	using simple_moves::SmoothSymmetricFragmentMover;
	using simple_moves::GunnCost;
	if ( option[ OptionKeys::abinitio::log_frags ].user() ) {
		if ( !option[ OptionKeys::abinitio::debug ] ) utility_exit_with_message( "apply option abinitio::log_frags always together with abinitio::debug!!!");
		bms = simple_moves::ClassicFragmentMoverOP( new simple_moves::LoggedFragmentMover( fragset_small, movemap ) );
		bml = simple_moves::ClassicFragmentMoverOP( new simple_moves::LoggedFragmentMover( fragset_large, movemap ) );
		sms = simple_moves::ClassicFragmentMoverOP( new SmoothFragmentMover( fragset_small, movemap, FragmentCostOP( new GunnCost ) ) );
	} else if ( option[ OptionKeys::abinitio::symmetry_residue ].user() ) {
		Size const sr (  option[ OptionKeys::abinitio::symmetry_residue ] );
		bms = simple_moves::ClassicFragmentMoverOP( new SymmetricFragmentMover( fragset_small, movemap, sr ) );
		bml = simple_moves::ClassicFragmentMoverOP( new SymmetricFragmentMover( fragset_large, movemap, sr ) );
		sms = simple_moves::ClassicFragmentMoverOP( new SmoothSymmetricFragmentMover( fragset_small, movemap, FragmentCostOP( new GunnCost ), sr ) );
	} else {
		bms = simple_moves::ClassicFragmentMoverOP( new ClassicFragmentMover( fragset_small, movemap ) );
		bml = simple_moves::ClassicFragmentMoverOP( new ClassicFragmentMover( fragset_large, movemap ) );
		sms = simple_moves::ClassicFragmentMoverOP( new SmoothFragmentMover ( fragset_small, movemap, FragmentCostOP( new GunnCost ) ) );
	}

	bms->set_end_bias( option[ OptionKeys::abinitio::end_bias ] ); //default is 30.0
	bml->set_end_bias( option[ OptionKeys::abinitio::end_bias ] );
	sms->set_end_bias( option[ OptionKeys::abinitio::end_bias ] );

	brute_move_small_ = bms;
	brute_move_large_ = bml;
	smooth_move_small_ = sms;

	using namespace core::pack::task;
	//init the packer
	pack_rotamers_ = minimization_packing::PackRotamersMoverOP( new protocols::minimization_packing::PackRotamersMover() );
	TaskFactoryOP main_task_factory( new TaskFactory );
	main_task_factory->push_back( operation::TaskOperationCOP( new operation::RestrictToRepacking ) );
	//main_task_factory->push_back( new operation::PreserveCBeta );
	pack_rotamers_->task_factory(main_task_factory);

	bSkipStage1_ = false;
	bSkipStage2_ = false;

	close_chbrk_ = false;

	stage4_cycles_pack_rate_ = 0.25;
}

/// @details Call parent's copy constructor and perform a shallow
/// copy of all the data.  NOTE: Shallow copy is only to preserve
/// behavior pre 9/7/2009 when the compiler-provided copy constructor
/// was being invoked.
ClassicAbinitio::ClassicAbinitio( ClassicAbinitio const & src ) :
	//utility::pointer::ReferenceCount(),
	Parent( src )
{
	stage1_cycles_ = src.stage1_cycles_;
	stage2_cycles_ = src.stage2_cycles_;
	stage3_cycles_ = src.stage3_cycles_;
	stage4_cycles_ = src.stage4_cycles_;
	stage5_cycles_ = src.stage5_cycles_;
	score_stage1_ = src.score_stage1_;
	score_stage2_ = src.score_stage2_;
	score_stage3a_ = src.score_stage3a_;
	score_stage3b_ = src.score_stage3b_;
	score_stage4_ = src.score_stage4_;
	score_stage4rot_ = src.score_stage4rot_;
	score_stage5_ = src.score_stage5_;
	apply_large_frags_ = src.apply_large_frags_;
	short_insert_region_ = src.short_insert_region_;
	just_smooth_cycles_ = src.just_smooth_cycles_;
	bQuickTest_ = src.bQuickTest_;
	close_chbrk_ = src.close_chbrk_;
	temperature_ = src.temperature_;
	movemap_ = src.movemap_;
	mc_ = src.mc_;
	brute_move_small_ = src.brute_move_small_;
	brute_move_large_ = src.brute_move_large_;
	smooth_move_small_ = src.smooth_move_small_;
	trial_large_ = src.trial_large_;
	trial_small_ = src.trial_small_;
	smooth_trial_small_ = src.smooth_trial_small_;
	total_trials_ = src.total_trials_;
	bSkipStage1_ = src.bSkipStage1_;
	bSkipStage2_ = src.bSkipStage2_;
	bSkipStage3_ = src.bSkipStage3_;
	bSkipStage4_ = src.bSkipStage4_;
	bSkipStage5_ = src.bSkipStage5_;
	recover_low_stages_ = src.recover_low_stages_;
}

/// @brief Explicit destructor is needed to destroy all the OPs
/// The compiler does all the work, but it requires that we place
/// the destructor in the .cc file.
ClassicAbinitio::~ClassicAbinitio() = default;

/// @brief setup moves, mc-object, scores
/// @details can't call this from constructor; virtual functions don't operate until construction has completed.

void
ClassicAbinitio::init( core::pose::Pose const& pose ) {
	// Parent::init( pose );
	set_defaults( pose );
	// bInitialized_ = true;
}

/// @brief ClassicAbinitio has virtual functions... use this to obtain a new instance
moves::MoverOP
ClassicAbinitio::clone() const
{
	return moves::MoverOP( new ClassicAbinitio( *this ) );
}


void ClassicAbinitio::apply( core::pose::Pose & pose ) {
	using namespace moves;
	using namespace scoring;
	using namespace scoring::constraints;

	Parent::apply( pose );
	if ( option[ OptionKeys::run::dry_run ]() ) return;
	//basic::prof_reset();
	total_trials_ = 0;
	
//.......................................................................................................
	lenth_of_sequence = pose.total_residue();
	
	Read_parameters();
	get_parameters(parametersMap);
	
	Read_distance();
	
	
	initPose = pose;
	vector<core::pose::Pose> storage_pose = Initial_population( pose );
	
	
	cout << "===========================Modal Exploration===========================" << endl;
	Multimodal_explore(storage_pose, G_explore*0.4);
	cout << "#######################################################################" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                     Modal Exploration Complete                      #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#######################################################################" << endl;
	cout << endl;

	
	cout << "===========================Modal Maintaining===========================" << endl;
	num_false = 0;
	bool can_forming;
	can_forming = Multimodal_maintain(storage_pose);
	
	while(can_forming == false){
		Multimodal_explore(storage_pose, G_explore*0.2);
		can_forming = Multimodal_maintain(storage_pose);
		
		num_false ++;
	}
	cout << "#######################################################################" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                     Modal Maintaining Complete                      #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#######################################################################" << endl;
	cout << endl;

	cout << "===========================Modal Exploitation==========================" << endl;
	Multimodal_enhance_c();
	cout << "#######################################################################" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                      Modal Exploitation Complete                    #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#                                                                     #" << endl;
	cout << "#######################################################################" << endl;
	cout << endl;
	final_enhence(Vec_modal_contact);
}

void 
ClassicAbinitio::Multimodal_explore(vector<core::pose::Pose> &storage_pose, int Gen){
///================================ @Multimodal_explore ================================
	KT_D = 1;
	frag_ = 9;
	FragAssem_ = brute_move_large_;
	

 	vector< vector<double> > DMscore_map( calculate_population_DMscore(storage_pose) );
	vector<double> storage_pose_Cscore( calculate_population_Dscore(storage_pose) );
	
	for(int g = 1; g <= Gen; g++){
		Multimodal_explore(storage_pose, storage_pose_Cscore, DMscore_map);
	}
}
	
bool 
ClassicAbinitio::Multimodal_maintain(vector<core::pose::Pose> &storage_pose){
///================================ @Multimodal_maintain ================================
	frag_ = 9;
	FragAssem_ = brute_move_large_;
	
	DMscore_and_cluster( storage_pose );
	
	double M_score = 0.5;
	vector< vector< vector<double> > > DMscore_map2;
	vector<core::pose::Pose> low_modal;
	vector<core::pose::Pose> top_modal;
	vector< vector<core::pose::Pose> > vec_top_modal;
	if(num_false <= 2){
		for(Size i = 0; i < best_every_cluster.size(); i++){
			vector<core::pose::Pose> population;
			for(Size j = 0; j < best_every_cluster[i].size(); j++){
				population.push_back( storage_pose[best_every_cluster[i][j]] );
			}
			DMscore_map2.push_back( calculate_population_DMscore(population) );
			
			double Num_DMscore = 0;
			for(Size j = 0; j < DMscore_map2[i].size(); j++)
				for(Size k = j+1; k < DMscore_map2[i][j].size(); k++)
					Num_DMscore += DMscore_map2[i][j][k];
			double ave_DMscore = Num_DMscore/(best_every_cluster[i].size()*(best_every_cluster[i].size()-1)*0.5);
			
			double Modal_score;
			if(best_every_cluster[i].size() >= 2)
				Modal_score = 2/(1/(log(best_every_cluster[i].size())/log(N)) + 1/ave_DMscore);
			else
				Modal_score = 0;
			
			if( Modal_score >= M_score ){
				vector<core::pose::Pose> ith_top_modal;
				for(Size j = 0; j < best_every_cluster[i].size(); j++){
					ith_top_modal.push_back( storage_pose[best_every_cluster[i][j]] );
					top_modal.push_back(storage_pose[best_every_cluster[i][j]]);
				}
				vec_top_modal.push_back(ith_top_modal); 
			}
			else{
				for(Size j = 0; j < best_every_cluster[i].size(); j++)
					low_modal.push_back( storage_pose[best_every_cluster[i][j]] );
			}
		}
	}
	else{
		vector< pair<double, int> > Modal_score_i;
		for(Size i = 0; i < best_every_cluster.size(); i++){
			vector<core::pose::Pose> population;
			for(Size j = 0; j < best_every_cluster[i].size(); j++){
				population.push_back( storage_pose[best_every_cluster[i][j]] );
			}
			DMscore_map2.push_back( calculate_population_DMscore(population) );
			
			double Num_DMscore = 0;
			for(Size j = 0; j < DMscore_map2[i].size(); j++)
				for(Size k = j+1; k < DMscore_map2[i][j].size(); k++)
					Num_DMscore += DMscore_map2[i][j][k];
			double ave_DMscore = Num_DMscore/(best_every_cluster[i].size()*(best_every_cluster[i].size()-1)*0.5);
			
			double Modal_score;
			if(best_every_cluster[i].size() >= 2)
				Modal_score = 2/(1/(log(best_every_cluster[i].size())/log(N)) + 1/ave_DMscore);
			else
				Modal_score = 0;
						
			Modal_score_i.push_back( make_pair(Modal_score, i) );
		}
		
		sort(Modal_score_i.begin(), Modal_score_i.end());
		for(Size m = best_every_cluster.size()-12; m < best_every_cluster.size(); m++ ){
			vector<core::pose::Pose> ith_top_modal;
			for(Size j = 0; j < best_every_cluster[Modal_score_i[m].second].size(); j++){
				ith_top_modal.push_back( storage_pose[best_every_cluster[Modal_score_i[m].second][j]] );
				top_modal.push_back(storage_pose[best_every_cluster[Modal_score_i[m].second][j]]);
			}
			
			vec_top_modal.push_back(ith_top_modal); 
		}
		for(Size m = 0; m < best_every_cluster.size()-12; m++ ){
			for(Size j = 0; j < best_every_cluster[Modal_score_i[m].second].size(); j++)
				low_modal.push_back( storage_pose[best_every_cluster[Modal_score_i[m].second][j]] );
		}
	}
	
	bool can_for(true);
	if(top_modal.size() >= storage_pose.size()/2 || num_false > 2 ){
		///@note =============================================Modal merging============================================	
		double Cut_DMscore = 0.5;
		Size gen_no_change = 0;
		for (int g = 1; g <= G_maintain*0.5; ++g){
			
			vector< vector<double> > true_modal_Cscore_mer;
			vector< vector<double> > true_modal_energy_mer;
			vector< vector< vector<double> > > DMscore_map_mer;
			///@note ######################### true-positive modal evolution #########################
			
			for(Size i = 0; i < vec_top_modal.size(); i++){
				vector<core::pose::Pose> population;
				for(Size j = 0; j < vec_top_modal[i].size(); j++)
					population.push_back( vec_top_modal[i][j] );
				
				true_modal_Cscore_mer.push_back( calculate_population_Dscore(population) );
				true_modal_energy_mer.push_back( calculate_population_energy(population) );
				
				Multimodal_maintain(population, true_modal_energy_mer[i], true_modal_Cscore_mer[i]);
				
				DMscore_map_mer.push_back( calculate_population_DMscore(population) );
			}
			//记录质心
			vector<core::pose::Pose> center_of_top_modal;
			for(Size i = 0; i < DMscore_map_mer.size(); i++){
				vector<double> sum_j; 
				for(Size j = 0; j < DMscore_map_mer[i].size(); j++){
					double sum_dmscore = 0;
					for(Size k = 0; k < DMscore_map_mer[i][j].size(); k++){
						sum_dmscore += DMscore_map_mer[i][j][k];
					}
					sum_j.push_back( sum_dmscore );
				}
				int max_sum_j_index = max_element(sum_j.begin(), sum_j.end()) - sum_j.begin();
				center_of_top_modal.push_back(vec_top_modal[i][max_sum_j_index]);
				
			}
			
			///@note ######################### false-positive modal evolution #########################
			Size last_gen_size = low_modal.size();
			vector<core::pose::Pose> merging_fail_pose;
			merging_fail_pose = Modal_merging(low_modal,top_modal, center_of_top_modal, vec_top_modal, Cut_DMscore);
			
			if(merging_fail_pose.size() == 0)
				break;
			
			vector<core::pose::Pose>().swap(low_modal);
			for(Size i = 0; i < merging_fail_pose.size(); i++)
				low_modal.push_back(merging_fail_pose[i]);
			
			
			if(low_modal.size() == last_gen_size)
				gen_no_change ++;
			else
				gen_no_change = 0;
			
			if(gen_no_change >= 10){
				Cut_DMscore -= 0.02;
				gen_no_change = 0;
			}
			
			if(g == G_maintain*0.5 && low_modal.size() != 0){
				for(Size i = 0; i < low_modal.size(); i++){
					vector<double> vec_trialPose;
					for(Size k = 0; k < center_of_top_modal.size(); k++)
						vec_trialPose.push_back( DM_score(low_modal[i], center_of_top_modal[k]) );
					
					int max_vec_trialPose_index = max_element(vec_trialPose.begin(), vec_trialPose.end())-vec_trialPose.begin();
					
					vec_top_modal[max_vec_trialPose_index].push_back(low_modal[i]);
					top_modal.push_back(low_modal[i]);
				}
			}
		}
		for(Size i = 0; i < top_modal.size(); i++)
			storage_pose[i] = top_modal[i];
		
		///@note ######################### all modal evolution #########################
		for (int g = 1; g <= G_maintain*0.5; ++g){
			vector< vector<double> > true_modal_Cscore_mer;
			vector< vector<double> > true_modal_energy_mer;
			for(Size i = 0; i < vec_top_modal.size(); i++){
				vector<core::pose::Pose> population;
				for(Size j = 0; j < vec_top_modal[i].size(); j++)
					population.push_back( vec_top_modal[i][j] );
				
				true_modal_Cscore_mer.push_back( calculate_population_Dscore(population) );
				true_modal_energy_mer.push_back( calculate_population_energy(population) );
				
				Multimodal_maintain(population, true_modal_energy_mer[i], true_modal_Cscore_mer[i]);
			}
		}
	}
	else{
		can_for = false;
	}
	
	Vec_modal_contact = vec_top_modal;
	
	return can_for;
}

void 
ClassicAbinitio::Multimodal_enhance_c(){
///================================ @Multimodal_enhance_c ================================
	KT_D = 1;
	frag_ = 3;
	FragAssem_ = brute_move_small_;
	
	for (Size m = 0; m < Vec_modal_contact.size(); m++){
		vector<core::pose::Pose> population;
		vector<double> vec_modal_energy;
		vector<double> vec_modal_Cscore;
		
		for(Size i = 0; i < Vec_modal_contact[m].size(); i++)
			population.push_back( Vec_modal_contact[m][i] );
		
		vec_modal_energy = calculate_population_energy(population);
		vec_modal_Cscore = calculate_population_Dscore(population);
		Ave_Dscore = Average_score_of_Populatin( vec_modal_Cscore );
		
		for (int g = 0; g < G_enhance_c; g++){
			Multimodal_enhance_c(population, vec_modal_energy, vec_modal_Cscore, g);
		}
		
		for(Size i = 0; i < Vec_modal_contact[m].size(); i++)
			Vec_modal_contact[m][i] = population[i];
	}
}

void 
ClassicAbinitio::Multimodal_explore(vector< core::pose::Pose >& storage_pose, vector< double >& population_Dscore, vector< vector<double> > &DMscore_map){
///=================================== @Multimodal_forming ================================================
	int NP_ = storage_pose.size();
	for ( int i = 0; i < NP_; ++i){
		int trialTimes ( 0 );
		core::pose::Pose trialPose;
	
		int base( numeric::random::rg().random_range(0, NP_ - 1) );
		core::pose::Pose base_pose( storage_pose[base] );
		
		///====================== @Fragment_recombination =====================
		int rand2( numeric::random::rg().random_range(0, NP_ - 1) );
		int rand3( numeric::random::rg().random_range(0, NP_ - 1) );
		while ( base == i )
			base = numeric::random::rg().random_range(0, NP_ - 1);
		while ( rand2 == base || rand2 == i )
			rand2 = numeric::random::rg().random_range(0, NP_ - 1);
		while ( rand3 == base || rand3 == rand2 || rand3 == i )
			rand3 = numeric::random::rg().random_range(0, NP_ - 1);
		
		core::pose::Pose rand2_pose( storage_pose[rand2] );
		core::pose::Pose rand3_pose( storage_pose[rand3] );
		
		int rand_w1( numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1) );
		int rand_w2( numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1) );
		int rand_w3( numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1) );
		while ( rand_w2 == rand_w1 )
			rand_w2 = numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1);
		while ( rand_w3 == rand_w1 || rand_w3 == rand_w2 )
			rand_w3 = numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1);
		
		for (int r = rand_w1; r < rand_w1 + frag_; ++r){
			base_pose.set_phi( r, rand2_pose.phi( r ) );
			base_pose.set_psi( r, rand2_pose.psi( r ) );
			base_pose.set_omega( r, rand2_pose.omega( r ) );
		}
		for (int r = rand_w2; r < rand_w2 + frag_; ++r){
			base_pose.set_phi( r, rand3_pose.phi( r ) );
			base_pose.set_psi( r, rand3_pose.psi( r ) );
			base_pose.set_omega( r, rand3_pose.omega( r ) );
		}
		for (int r = rand_w3; r < rand_w3 + frag_; ++r){
			base_pose.set_phi( r, storage_pose[i].phi( r ) );
			base_pose.set_psi( r, storage_pose[i].psi( r ) );
			base_pose.set_omega( r, storage_pose[i].omega( r ) );
		}
	
	
		///==================== @Fragment_assemble ========================
		double baseEnergy( (*score_stage4_)(base_pose) );
		double trialEnergy;
		while ( trialTimes < 200 ){
			trialPose = base_pose;
			FragAssem_->apply( trialPose );
			trialEnergy = (*score_stage4_)(trialPose);
			if ( boltzmann_accept(baseEnergy, trialEnergy) )
				break;
			++trialTimes;
		}
		if ( trialTimes >= 200 ){
			trialPose = base_pose;
			trialEnergy = baseEnergy;
		}
		
		///================= @Individual_Crowding ================ @Individual_Crowding ================
		double trialCscore;
		trialCscore = Distance_score(trialPose);
		
		//============寻找vec_trialPose和DMscore_map中DMscore最大的============================
		vector<double> vec_trialPose;
		for(int k = 0; k < NP_; k++){
			vec_trialPose.push_back( DM_score(trialPose, storage_pose[k]) );
		}
		

		vector< pair<double, int> > CMscore_and_num;
		for( int i = 0; i < NP_; i++ )
			CMscore_and_num.push_back( make_pair(vec_trialPose[i], i) );
		
		sort(CMscore_and_num.begin(), CMscore_and_num.end());
		
		vector<int> rand_integer;
		for(int i = 0; i < Crowding_radius_exp; i++){
			rand_integer.push_back( i );
		}
		
		random_shuffle( rand_integer.begin(), rand_integer.end() );
		
		for(int r = 0; r < Crowding_radius_exp; r++){
			int crow_index = CMscore_and_num[NP_-1-rand_integer[r]].second;
			if ( boltzmann_accept(population_Dscore[crow_index], trialCscore, KT_D) ){
				for(int k = 0; k < NP_; k++){
					DMscore_map[crow_index][k] = vec_trialPose[k];
					DMscore_map[k][crow_index] = vec_trialPose[k];
				}
				DMscore_map[crow_index][crow_index] = 1;
				
				storage_pose[crow_index] = trialPose;
				population_Dscore[crow_index] = trialCscore;
				
				break;
			}
		}
		
	}
}

void 
ClassicAbinitio::Multimodal_maintain(vector< core::pose::Pose >& population, vector< double >& population_energy, vector< double >& population_Dscore){
	KT_E = 1;
	KT_D = 1;
	int NP_ = population.size();
	for ( int i = 0; i < NP_; ++i){
		int base = 0;
		core::pose::Pose base_pose;
		
		base = numeric::random::rg().random_range(0, NP_ - 1);
		base_pose = population[base];
		
		if( NP_ >= 4 ){
			///====================== @Fragment_recombination =====================
			int rand2( numeric::random::rg().random_range(0, NP_ - 1) );
			int rand3( numeric::random::rg().random_range(0, NP_ - 1) );
			while ( base == i )
				base = numeric::random::rg().random_range(0, NP_ - 1);
			while ( rand2 == base || rand2 == i )
				rand2 = numeric::random::rg().random_range(0, NP_ - 1);
			while ( rand3 == base || rand3 == rand2 || rand3 == i )
				rand3 = numeric::random::rg().random_range(0, NP_ - 1);
			
			core::pose::Pose rand2_pose( population[rand2] );
			core::pose::Pose rand3_pose( population[rand3] );
			
			int rand_w1( numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1) );
			int rand_w2( numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1) );
			int rand_w3( numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1) );
			while ( rand_w2 == rand_w1 )
				rand_w2 = numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1);
			while ( rand_w3 == rand_w1 || rand_w3 == rand_w2 )
				rand_w3 = numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1);
			
			for (int r = rand_w1; r < rand_w1 + frag_; ++r){
				base_pose.set_phi( r, rand2_pose.phi( r ) );
				base_pose.set_psi( r, rand2_pose.psi( r ) );
				base_pose.set_omega( r, rand2_pose.omega( r ) );
			}
			for (int r = rand_w2; r < rand_w2 + frag_; ++r){
				base_pose.set_phi( r, rand3_pose.phi( r ) );
				base_pose.set_psi( r, rand3_pose.psi( r ) );
				base_pose.set_omega( r, rand3_pose.omega( r ) );
			}
			for (int r = rand_w3; r < rand_w3 + frag_; ++r){
				base_pose.set_phi( r, population[i].phi( r ) );
				base_pose.set_psi( r, population[i].psi( r ) );
				base_pose.set_omega( r, population[i].omega( r ) );
			}
		}
		
		///==================== @Fragment_assemble ========================
		double baseEnergy( (*score_stage4_)(base_pose) );
		core::pose::Pose trialPose;
		double trialEnergy;
		int trialTimes ( 0 );
		while ( trialTimes < 200 ){
			trialPose = base_pose;
			FragAssem_->apply( trialPose );
			trialEnergy = (*score_stage4_)(trialPose);
			if ( boltzmann_accept(baseEnergy, trialEnergy, KT_E) )
				break;
			++trialTimes;
		}
		if ( trialTimes >= 200 ){
			trialPose = base_pose;
			trialEnergy = baseEnergy;
		}
		
		///================= @Individual_Crowding ================ @Individual_Crowding ================
		double trialenergy( (*score_stage4_)(trialPose) );
		double trialCscore( Distance_score(trialPose) );
		
		vector<double> vec_trialPose;
		for(int k = 0; k < NP_; k++)
			vec_trialPose.push_back( DM_score(trialPose, population[k]) );
		
		
		vector< pair<double, int> > CMscore_and_num;
		for( int i = 0; i < NP_; i++ )
			CMscore_and_num.push_back( make_pair(vec_trialPose[i], i) );
		sort(CMscore_and_num.begin(), CMscore_and_num.end());
		
		int Crowding_radius = Crowding_radius_mer;
		if(NP_ < Crowding_radius_mer){
			Crowding_radius = NP_;
		}
		
		vector<int> rand_integer;
		for(int i = 0; i < Crowding_radius; i++){
			rand_integer.push_back( i );
		}
		
		random_shuffle( rand_integer.begin(), rand_integer.end() );
		
		for(int r = 0; r < Crowding_radius; r++){
			int crow_index = CMscore_and_num[NP_-1-rand_integer[r]].second;
			if ( boltzmann_accept(population_energy[crow_index], trialenergy, KT_E) ){
				population[crow_index] = trialPose;
				population_energy[crow_index] = trialEnergy;
				population_Dscore[crow_index] = trialCscore;
				
			}
			else if ( boltzmann_accept(population_Dscore[crow_index], trialCscore, KT_D) ){
				population[crow_index] = trialPose;
				population_energy[crow_index] = trialEnergy;
				population_Dscore[crow_index] = trialCscore;
				
			}
		}
	}
}

vector<core::pose::Pose> 
ClassicAbinitio::Modal_merging(vector< core::pose::Pose >& population, vector< core::pose::Pose >& top_population, 
			       vector<core::pose::Pose> &center_population, vector< vector<core::pose::Pose> > &vec_top_modal, double &cut_dmscore){
	KT_E = 5;
	int NP_ = population.size();
	int NP_top = top_population.size();
	vector<core::pose::Pose> merging_fail;
	for ( int i = 0; i < NP_; ++i){
		int base = 0;
		core::pose::Pose base_pose;
		base_pose = population[i];
		
		///====================== @Fragment_recombination =====================
		int rand2( numeric::random::rg().random_range(0, NP_top - 1) );
		int rand3( numeric::random::rg().random_range(0, NP_top - 1) );
		while ( base == i )
			base = numeric::random::rg().random_range(0, NP_top - 1);
		while ( rand2 == base || rand2 == i )
			rand2 = numeric::random::rg().random_range(0, NP_top - 1);
		while ( rand3 == base || rand3 == rand2 || rand3 == i )
			rand3 = numeric::random::rg().random_range(0, NP_top - 1);
		
		core::pose::Pose rand2_pose( top_population[rand2] );
		core::pose::Pose rand3_pose( top_population[rand3] );
		
		int rand_w1( numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1) );
		int rand_w2( numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1) );
		while ( rand_w2 == rand_w1 )
			rand_w2 = numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1);
		
		for (int r = rand_w1; r < rand_w1 + frag_; ++r){
			base_pose.set_phi( r, rand2_pose.phi( r ) );
			base_pose.set_psi( r, rand2_pose.psi( r ) );
			base_pose.set_omega( r, rand2_pose.omega( r ) );
		}
		for (int r = rand_w2; r < rand_w2 + frag_; ++r){
			base_pose.set_phi( r, rand3_pose.phi( r ) );
			base_pose.set_psi( r, rand3_pose.psi( r ) );
			base_pose.set_omega( r, rand3_pose.omega( r ) );
		}
		
		///==================== @Fragment_assemble ========================
		double baseEnergy( (*score_stage4_)(base_pose) );
		core::pose::Pose trialPose;
		double trialEnergy;
		int trialTimes ( 0 );
		while ( trialTimes < 200 ){
			trialPose = base_pose;
			FragAssem_->apply( trialPose );
			trialEnergy = (*score_stage4_)(trialPose);
			if ( boltzmann_accept(baseEnergy, trialEnergy, KT_E) )
				break;
			++trialTimes;
		}
		if ( trialTimes >= 200 ){
			trialPose = base_pose;
			trialEnergy = baseEnergy;
		}
		
		vector<double> vec_trialPose;
		for(Size k = 0; k < center_population.size(); k++)
			vec_trialPose.push_back( DM_score(trialPose, center_population[k]) );
		double max_vec_trialPose = *max_element(vec_trialPose.begin(), vec_trialPose.end());
		int max_vec_trialPose_index = max_element(vec_trialPose.begin(), vec_trialPose.end())-vec_trialPose.begin();
		
		if(max_vec_trialPose >= cut_dmscore){
			vec_top_modal[max_vec_trialPose_index].push_back(trialPose);
			top_population.push_back(trialPose);
		}
		else
			merging_fail.push_back(population[i]);
	}
	
	
	return merging_fail;
}

void 
ClassicAbinitio::Multimodal_enhance_c(vector<core::pose::Pose> &population, vector<double> &population_energy, vector<double> &population_Dscore, int enhance_g){
	int NP_ = population.size();
	for ( int i = 0; i < NP_; ++i){
		///====================== @Fragment_recombination =====================
		int base( numeric::random::rg().random_range(0, NP_ - 1) );
		while ( base == i )
			base = numeric::random::rg().random_range(0, NP_ - 1);
		core::pose::Pose base_pose( population[base] );
		if( NP_ >= 4 ){
			///====================== @Fragment_recombination =====================
			int rand2( numeric::random::rg().random_range(0, NP_ - 1) );
			int rand3( numeric::random::rg().random_range(0, NP_ - 1) );
			while ( rand2 == base || rand2 == i )
				rand2 = numeric::random::rg().random_range(0, NP_ - 1);
			while ( rand3 == base || rand3 == rand2 || rand3 == i )
				rand3 = numeric::random::rg().random_range(0, NP_ - 1);
			
			core::pose::Pose rand2_pose( population[rand2] );
			core::pose::Pose rand3_pose( population[rand3] );
			
			int rand_w1( numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1) );
			int rand_w2( numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1) );
			int rand_w3( numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1) );
			while ( rand_w2 == rand_w1 )
				rand_w2 = numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1);
			while ( rand_w3 == rand_w1 || rand_w3 == rand_w2 )
				rand_w3 = numeric::random::rg().random_range(1, lenth_of_sequence - frag_ + 1);
			
			for (int r = rand_w1; r < rand_w1 + frag_; ++r){
				base_pose.set_phi( r, rand2_pose.phi( r ) );
				base_pose.set_psi( r, rand2_pose.psi( r ) );
				base_pose.set_omega( r, rand2_pose.omega( r ) );
			}
			for (int r = rand_w2; r < rand_w2 + frag_; ++r){
				base_pose.set_phi( r, rand3_pose.phi( r ) );
				base_pose.set_psi( r, rand3_pose.psi( r ) );
				base_pose.set_omega( r, rand3_pose.omega( r ) );
			}
			for (int r = rand_w3; r < rand_w3 + frag_; ++r){
				base_pose.set_phi( r, population[i].phi( r ) );
				base_pose.set_psi( r, population[i].psi( r ) );
				base_pose.set_omega( r, population[i].omega( r ) );
			}
		}
		
		///==================== @Fragment_assemble ========================
		double baseEnergy( (*score_stage4_)(base_pose) );
		core::pose::Pose trialPose;
		double trialEnergy;
		int trialTimes ( 0 );
		while ( trialTimes < 200 ){
			trialPose = base_pose;
			FragAssem_->apply( trialPose );
			trialEnergy = (*score_stage4_)(trialPose);
			if ( boltzmann_accept(baseEnergy, trialEnergy) )
				break;
			++trialTimes;
		}
		if ( trialTimes >= 200 ){
		//	cout << "### Fragment Assemble Failed Afeter 200 trials!!!" << endl;
			trialPose = base_pose;
			trialEnergy = baseEnergy;
		}
		///================= @Individual_Selection ================ @Individual_Selection ================
		double trialCscore( Distance_score(trialPose) );
		TrialEnergy_ = trialEnergy;
		
		if(enhance_g <= (G_enhance_c/8)*5){
			if ( find_if( population_energy.begin(), population_energy.end(), isSimilar ) == population_energy.end() ){
				if ( boltzmann_accept(population_Dscore[i], trialCscore, KT_D) ){
					population[i] = trialPose;
					population_energy[i] = trialEnergy;
					Ave_Dscore = Average_score_of_Populatin( population_Dscore[i], trialCscore, Ave_Dscore, NP_ );
					population_Dscore[i] = trialCscore;
				}
				else if ( boltzmann_accept(Ave_Dscore, trialCscore, KT_D) && boltzmann_accept(population_energy[i], trialEnergy) ){
					population[i] = trialPose;
					population_energy[i] = trialEnergy;
					Ave_Dscore = Average_score_of_Populatin( population_Dscore[i], trialCscore, Ave_Dscore, NP_ );
					population_Dscore[i] = trialCscore;
				}
			}
		}
		else {
			int max_Cscore_Pose_index = max_element(population_Dscore.begin(), population_Dscore.end())-population_Dscore.begin();
			if ( find_if( population_energy.begin(), population_energy.end(), isSimilar ) == population_energy.end() ){
				if ( boltzmann_accept(population_Dscore[max_Cscore_Pose_index], trialCscore, KT_D) ){
					population[max_Cscore_Pose_index] = trialPose;
					population_energy[max_Cscore_Pose_index] = trialEnergy;
					population_Dscore[max_Cscore_Pose_index] = trialCscore;
				}
				else if ( boltzmann_accept(population_energy[i], trialEnergy) ){
					population[i] = trialPose;
					population_energy[i] = trialEnergy;
					population_Dscore[i] = trialCscore;
				}
			}
		}
	}
}

void 
ClassicAbinitio::Kmediods(int K, vector< vector<double> > Distance_Matrix){
	if(K >= N ){
		cout << "=========================The population is too small !!!" << endl;
		exit(0);
	}
	
	vector<int> mediods;
	vector< vector<int> > every_cluster;
	///@note the number of clustring points.
	int NN = Distance_Matrix.size();
	///@note store K clustring centers.
	mediods.resize(K,0);
	///@brief randomly generates the cluster centers.
	for ( Size i = 0; i < mediods.size(); ++i ){
		int rand_mediod( 0 );
		while ( 1 ){
			rand_mediod = numeric::random::rg().random_range(0, NN-1);
			Size j = 0;
			for (; j < mediods.size(); ++j){
				if ( rand_mediod == mediods[j] )
					break;
			}
			if ( j == mediods.size() ){
				mediods[i] = rand_mediod;
				break;
			}
		}
	}
	
	///@note  store the cluster to which the point belong for each clustring points.
	vector< unsigned int > clusterAssement( NN, 0 );
	///@brief cluster n points according to the cluster centers.
	for ( int i = 0; i < NN; ++i ){
		///@note store closest cluster and the distance.
		pair<int, double> cluster_distance( make_pair(0, 0) );
		for ( Size k = 0; k < mediods.size(); ++k ){
			if( mediods[k]==i ){
				cluster_distance.first  = k;
				cluster_distance.second = Distance_Matrix[i][mediods[k]];
				break;
			}
			else if( Distance_Matrix[i][mediods[k]] < cluster_distance.second || k == 0 ){
				cluster_distance.first  = k;
				cluster_distance.second = Distance_Matrix[i][mediods[k]];
			}
		}
		clusterAssement[i] = cluster_distance.first;
	}
	
	///@brief avoiding have cluster with no points.
	int  num1 = 0;
	while ( num1 < 50 ){
		num1++;
		unsigned int k( 0 );
		for ( ; k < mediods.size(); ++k ){
			///@note extract all points belongs to cluster k.
			int count_points = 0;
			vector<int> points_in_cluster;
			for ( int i = 0; i < NN; ++i ){
				if ( clusterAssement[i] == k ){
					points_in_cluster.push_back( i );
					++count_points;
				}
			}
			if ( count_points <= 1 ){
				int rand_mediod( 0 );
				while ( 1 ){
					rand_mediod = numeric::random::rg().random_range(0, NN-1);
					Size j( 0 );
					for (; j < mediods.size(); ++j){
						if ( rand_mediod == mediods[j] )
							break;
					}
					if ( j == mediods.size() ){
						mediods[k] = rand_mediod;
						break;
					}
				}
				break;
			}
		}
		if ( k == mediods.size() )
			break;
		else{
			///@brief cluster n points according to the cluster centers.
			for ( int i = 0; i < NN; ++i ){
				///@note store closest cluster and the distance.
				pair<int, double> cluster_distance( make_pair(0, 0) );
				for ( Size k = 0; k < mediods.size(); ++k ){
				  
					if( mediods[k]==i ){
						cluster_distance.first  = k;
						cluster_distance.second = Distance_Matrix[i][mediods[k]];
						break;
					}
					else if ( Distance_Matrix[i][mediods[k]] < cluster_distance.second || k == 0 ){
						cluster_distance.first  = k;
						cluster_distance.second = Distance_Matrix[i][mediods[k]];
					}
				}
				clusterAssement[i] = cluster_distance.first;
			}
		}
		//eg: 聚类中心为2,3,8,clusterAssement[0]=3,clusterAssement[1]=2,clusterAssement[2]=3,clusterAssement[3]=8,...........
	}
	if ( num1 >=50 ){
		cout << "have cluster with 0 point, reselect a new cluster center!!!" << endl;
		exit(0);
	}
	
	///@note 记录类心是否发生改变,以及类心更新次数
	int  num2 = 0;
	///@brief iteration updata
	while ( num2 < 50){
		num2++;

		///@brief update cluster centers
// 		pair<int, int> max_cluster( make_pair(0, 0) );
		///max_cluster 记录最大的类和最大类中的个体数
		vector<int> count_points;
		count_points.resize(mediods.size(),0);
		for ( unsigned int k = 0; k < mediods.size(); ++k ){
			///@note extract all points belongs to cluster k.
			vector<int> points_in_cluster;
			for ( int i = 0; i < NN; ++i ){
				if ( clusterAssement[i] == k ){
					points_in_cluster.push_back( i );
					++count_points[k];
				}
			}//points_in_cluster 中存放这第k簇中的个体的下标
			
			///@note store cluster center and the total distance to all points.
			pair<int, double> center_distance(0, 0);
			///@note update
			for ( Size m = 0; m < points_in_cluster.size(); ++m ){
				double total_distance( 0 );
				for ( Size n = 0; n < points_in_cluster.size(); ++n )
					total_distance += Distance_Matrix[points_in_cluster[m]][points_in_cluster[n]];
				if ( total_distance < center_distance.second || m == 0 ){
					center_distance.first = points_in_cluster[m];
					center_distance.second = total_distance;
				}
			}
			mediods[k] = center_distance.first;
		}
		
		int count_of_state_transition = 0;
		///@brief update cluster (reclassification).
		for ( int i = 0; i < NN; ++i ){
			///@note store closest cluster and the distance.
			pair<unsigned int, double> cluster_distance( make_pair(0, 0) );
			for ( Size k = 0; k < mediods.size(); ++k ){
				if( mediods[k]==i ){
					cluster_distance.first  = k;
					cluster_distance.second = Distance_Matrix[i][mediods[k]];
					break;
				}
				else if ( Distance_Matrix[i][mediods[k]] < cluster_distance.second || k == 0 ){
					cluster_distance.first  = k;
					cluster_distance.second = Distance_Matrix[i][mediods[k]];
				}
			}
			if ( cluster_distance.first != clusterAssement[i] ){
				++count_of_state_transition;
				clusterAssement[i] = cluster_distance.first;
			}
		}		
		///@note if cluster is not change, break.
		if ( count_of_state_transition == 0 ){
			break;
		}
	}
	

	every_cluster.clear();
	every_cluster.resize( mediods.size() );
	for(int i=0; i<NN; i++){
		every_cluster[clusterAssement[i]].push_back(i);
	}
	
	
	int more_than_one = 0;
	double mean_num_distance = 0;
	for(Size i = 0; i < mediods.size(); i++){
		vector<double> vec_every_cluster;
		for(Size j = 0; j < every_cluster[i].size(); j++){
			vec_every_cluster.push_back( every_cluster[i][j] );
		}
		if( vec_every_cluster.size() > 1 ){
			more_than_one ++;
			double num_distance = 0;
			int num = 0;
			for( unsigned int k = 0; k < vec_every_cluster.size()-1; k++){
				for( unsigned int f = k+1; f < vec_every_cluster.size(); f++){
					num_distance += Distance_Matrix[vec_every_cluster[k]][vec_every_cluster[f]];
					num ++;
				}
			}
			mean_num_distance += num_distance/num;
		}
	}
	if( mean_num_distance/more_than_one < temp_mean_num_distance || K_Kmediods == 0 ){
		temp_mean_num_distance = mean_num_distance/more_than_one;
		best_every_cluster = every_cluster;
	}
}

void 
ClassicAbinitio::Map_matrix(vector< vector<double> > Map_matrix_pose){
	vector< vector<double> >temp_map_matrix_pose_100( N,vector<double>(N,0) );
	temp_map_matrix_pose_100 = Map_matrix_pose;
	
	double temp=0;
	for(int m=0; m<N; m++){
		for(int i=0; i<N; i++){
			for(int k=0; k<N-i-1; k++){
				if(Map_matrix_pose[m][k] > Map_matrix_pose[m][k+1]){
					temp = Map_matrix_pose[m][k];
					Map_matrix_pose[m][k] = Map_matrix_pose[m][k+1];
					Map_matrix_pose[m][k+1] = temp;
				}
			}
		}
	}
	
//////------------------------------------------------------方差----------------------------------------------------------
	//m阶距的和
	vector<double> num(N-1,0);
	for(int m=0; m<N-1; m++){
		for(int k=0; k<N; k++){
			num[m] = num[m]+Map_matrix_pose[k][m+1];
		}
	}
	//m阶距的平方和
	vector<double> mean_square(N-1,0);
	for(int m=0; m<N-1; m++){
		for(int k=0; k<N; k++){
			mean_square[m] = mean_square[m]+pow(Map_matrix_pose[k][m+1],2);
		}
	}
	//方差（平方均值-均值的平方）
	vector<double> temp_variance(N-1,0);
	vector<double> variance(N-1,0);
	for(int m=0; m<N-1; m++){
		for(int k=0; k<=m; k++){
			temp_variance[m] = temp_variance[m]+( mean_square[k]/N - pow(num[k]/N,2) );
		}
		variance[m] = temp_variance[m]/(m+1);
	}
	for(int m = 0; m < N-1; m++){
		variance[m] = 10000*(temp_variance[m]/(m+1));
	}
////------------------------------------------------------确定K个类----------------------------------------------------------
	double w = 3;
	int KK = 1;
	vector<double> variance_cha(N-2,0);
	vector<double> variance_dif(N-2,0);

	for(int m=0; m<N-2; m++){
		variance_cha[m] = fabs(variance[m+1]-variance[m]);
	}

	for(int m=0; m<N-3; m++){
		variance_dif[m] = fabs( fabs(variance[m+2]-variance[m+1]) - fabs(variance[m+1]-variance[m]) );
	}
	
//差的差归一化处理
// 	int Nor_KK = 0;
// 	double max_variance_dif = *max_element(variance_dif.begin(), variance_dif.end());
// 	double min_variance_dif = *min_element(variance_dif.begin(), variance_dif.end());
// 	for(int m=0; m<N-3; m++){
// 		double normalization_dif = (variance_dif[m]-min_variance_dif)/(max_variance_dif-min_variance_dif);
// 		if(normalization_dif > 0.2)
// 			Nor_KK++;
// 	}
	
	for(int m=0; m<N-4; m++){
		if(variance_dif[m+1]-variance_dif[m] > w * variance_dif[m]){
			KK = KK + 1;
		}
	}
	while(KK < 20){
		w -=0.2;
		for(int m=0; m<N-4; m++){
			if(variance_dif[m+1]-variance_dif[m] > w * variance_dif[m]){
				KK = KK + 1;
			}
		}
	}
	
	K_Kmediods = 0;
	temp_mean_num_distance = 0;
	for(int i = 0; i < 10; i++){
		Kmediods( KK, temp_map_matrix_pose_100 );
		K_Kmediods ++;
	}
}

void 
ClassicAbinitio::DMscore_and_cluster(vector<core::pose::Pose> &storage_pose){
	vector< vector<double> >TM_similar_score_map( storage_pose.size(),vector<double>(storage_pose.size(),0) );
	for(Size k = 0; k < storage_pose.size(); k++){
		for(Size f = k; f < storage_pose.size(); f++){
			double TM_two_similar_score;
			TM_two_similar_score = DM_score( storage_pose[k], storage_pose[f] );
			TM_similar_score_map[k][f] = TM_two_similar_score;
			TM_similar_score_map[f][k] = TM_two_similar_score;
		}
	}
	
	vector< vector<double> >map_matrix_pose_100( storage_pose.size(),vector<double>(storage_pose.size(),0) );
	for(Size m=0; m<storage_pose.size(); m++){
		for(Size k=0; k<storage_pose.size(); k++){
			map_matrix_pose_100[m][k] = 1 - TM_similar_score_map[m][k];
		}
	}
	
	Map_matrix( map_matrix_pose_100 );
}


// double 
// ClassicAbinitio::DM_score( core::pose::Pose &pose, core::pose::Pose &tempPose ){
// 	double score_sum = 0;
// 	double dm_score = 0;
// 	double num = 0;
// 	double d0 = 0;
// 	float varepsilon = 0.001;
// 	for(int i = 0; i < lenth_of_sequence; i++){
// 		for(int j = i+1; j < lenth_of_sequence; j++){
// 			Vector a1 = (pose.residue(i+1).atom("CA").xyz());
// 			Vector a2 = (pose.residue(j+1).atom("CA").xyz());
// 			
// 			Vector b1 = (tempPose.residue(i+1).atom("CA").xyz());
// 			Vector b2 = (tempPose.residue(j+1).atom("CA").xyz());
// 			
// 			d0 = log(varepsilon + fabs(i-j));
// 			double di = fabs( a1.distance(a2) - b1.distance(b2) );
// 			if( fabs(i-j) >= 3 ){
// 				score_sum += 1/(1 + pow(di/d0, 2));
// 				num += 1;
// 			}
// 		}
// 	}
// 	dm_score = score_sum/num;
// 	
// 	return dm_score;
// }

//快速版评估模型：相隔一个残基计算一次，不影响评估准确性以及预测精度． //推荐使用
double 
ClassicAbinitio::DM_score( core::pose::Pose &pose, core::pose::Pose &tempPose ){
	double score_sum = 0;
	double dm_score = 0;
	double num = 0;
	double d0 = 0;
	float varepsilon = 0.001;
	for(int i = 0; i < lenth_of_sequence;){
		for(int j = i+2; j < lenth_of_sequence;){
			Vector a1 = (pose.residue(i+1).atom("CA").xyz());
			Vector a2 = (pose.residue(j+1).atom("CA").xyz());
			
			Vector b1 = (tempPose.residue(i+1).atom("CA").xyz());
			Vector b2 = (tempPose.residue(j+1).atom("CA").xyz());
			
			d0 = log(varepsilon + fabs(i-j));
			double di = fabs( a1.distance(a2) - b1.distance(b2) );
			if( fabs(i-j) >= 3 ){
				score_sum += 1/(1 + pow(di/d0, 2));
				num += 1;
			}
			j = j+2;
		}
		i = i+2;
	}
	dm_score = score_sum/num;
	
	return dm_score;
}

bool 
ClassicAbinitio::boltzmann_accept(double const& targetEnergy, double const& trialEnergy){
	if ( trialEnergy <= targetEnergy )
		return true;
	else{
		if ( exp( -(trialEnergy - targetEnergy) * 0.5 ) < numeric::random::rg().uniform() )
			return false;
		else
			return true;
	}
}

bool 
ClassicAbinitio::boltzmann_accept(const double& targetEnergy, const double& trialEnergy, const double& recipocal_KT){
	double KT_ = 1/recipocal_KT;
	if ( trialEnergy <= targetEnergy )
		return true;
	else{
		if ( exp( -(trialEnergy - targetEnergy) * KT_ ) < numeric::random::rg().uniform() )
			return false;
		else
			return true;
	}
}

vector< double > 
ClassicAbinitio::calculate_population_energy(vector< core::pose::Pose >& population){
	vector<double> population_energy;
	
	for (Size i = 0; i < population.size(); ++i)
		population_energy.push_back( (*score_stage4_)(population[i]) );
	
	if ( population_energy.size() != population.size() ){
		cout << "ERROR!!!\n" << "ERROR : calculate population energy failed!!!\n" << "ERROR!!!" << endl;
		exit(0);
	}
	return population_energy;
}

vector< double > 
ClassicAbinitio::calculate_population_Dscore(vector< core::pose::Pose >& population){
	vector<double> population_Dscore;
	
	for (Size i = 0; i < population.size(); ++i)
		population_Dscore.push_back( Distance_score( population[i] ) );
		
	if ( population_Dscore.size() != population.size() ){
		cout << "ERROR!!!\n" << "ERROR : calculate population Cscore failed!!!\n" << "ERROR!!!" << endl;
		exit(0);
	}
	return population_Dscore;
}

vector< vector<double> >
ClassicAbinitio::calculate_population_DMscore(vector<core::pose::Pose> &storage_pose){
	int N_ = storage_pose.size();
	vector< vector<double> >DMscore_map( N_,vector<double>(N_,0) );
	for(int k = 0; k < N_; k++){
		for(int f = k; f < N_; f++){
			double TM_two_similar_score;
			TM_two_similar_score = DM_score( storage_pose[k], storage_pose[f] );
			DMscore_map[k][f] = TM_two_similar_score;
			DMscore_map[f][k] = TM_two_similar_score;
		}
	}
	return DMscore_map;
}

double 
ClassicAbinitio::Average_score_of_Populatin(vector< double >& population_score){
	double total_score( 0 );
	double NP_reciprocal = (double)1 / population_score.size();
	for (Size i = 0; i < population_score.size(); ++i)
		total_score += population_score[i];
	double average_score( total_score * NP_reciprocal );
	
	return average_score;
}

double 
ClassicAbinitio::Average_score_of_Populatin(double& old_score, double& new_score, double& average_score, const Size& population_size){
	double NP_reciprocal = (double)1 / population_size;
	return average_score + (new_score - old_score) * NP_reciprocal;
}

void 
ClassicAbinitio::Read_distance(){
	ifstream Rc("./input_files/distance");
	string distance_line;
	int num = 0;
	
	while ( getline(Rc, distance_line) ){
		
		istringstream Data(distance_line);
		DistanceType distance;
		
		Data >> distance.re1 >> distance.re2 >> distance.dist >> distance.confidence;
		vec_distance.push_back(distance);
	}
}

double 
ClassicAbinitio::Distance_score(core::pose::Pose& pose){
	double Dis_score = 0;
	double confidence_all = 0;
	double d0 = 0;
	double di = 0;
	for(Size k = 0; k < vec_distance.size(); k++){
		double C_distance = 0;
		if ( pose.residue(vec_distance[k].re1).name3() == "GLY" || pose.residue(vec_distance[k].re2).name3() == "GLY" ){
			numeric::xyzVector<Real> CA1 = pose.residue(vec_distance[k].re1).xyz("CA");
			numeric::xyzVector<Real> CA2 = pose.residue(vec_distance[k].re2).xyz("CA");
			C_distance = CA1.distance(CA2);
		}else{
			numeric::xyzVector<Real> CB1 = pose.residue(vec_distance[k].re1).xyz("CB");
			numeric::xyzVector<Real> CB2 = pose.residue(vec_distance[k].re2).xyz("CB");
			C_distance = CB1.distance(CB2);
		}
		if( fabs(vec_distance[k].re1 - vec_distance[k].re2) >= 3 ){
			d0 = log(fabs(vec_distance[k].re1 - vec_distance[k].re2));
			di = fabs(C_distance - vec_distance[k].dist);
			Dis_score += vec_distance[k].confidence*(di/d0);
		}
	}
	return Dis_score;
}


void 
ClassicAbinitio::final_enhence(vector< vector<core::pose::Pose> > Vec_modal){
	//对类进行排序
	vector< pair<double, int> > modal_mean_contact;
	for(Size i = 0; i < Vec_modal.size(); i++){
		double distance_num = 0;
		for(Size j = 0; j < Vec_modal[i].size(); j++){
			distance_num += Distance_score(Vec_modal[i][j]);
		}
		modal_mean_contact.push_back( make_pair(distance_num/Vec_modal[i].size(), i) );
	}
	
	sort(modal_mean_contact.begin(), modal_mean_contact.end());
	
	//选择类中能量最低的
	vector<core::pose::Pose> sort_score_pose;
	for(Size i = 0; i < modal_mean_contact.size(); i++){
		vector<double> last_modal_contact;
		for(Size j = 0; j < Vec_modal[modal_mean_contact[i].second].size(); j++)
			last_modal_contact.push_back( Distance_score(Vec_modal[modal_mean_contact[i].second][j]) );
		
		int min_contact_index = min_element(last_modal_contact.begin(), last_modal_contact.end()) - last_modal_contact.begin();
		sort_score_pose.push_back(Vec_modal[modal_mean_contact[i].second][min_contact_index]);
	}
	
	for(Size i = 0; i < 5; i++){
		Enhance_population.push_back(sort_score_pose[i]);
	}
	
	Enhance();
	
	vector< pair<double, int> > vec_dscore_i;
	for(Size i = 0; i < Enhance_population.size(); i++){
		vec_dscore_i.push_back( make_pair( Distance_score(Enhance_population[i]), i) );
		sort(vec_dscore_i.begin(), vec_dscore_i.end());
	}
	
	for(Size i = 0; i < Enhance_population.size(); i++){
		stringstream ss;
		string str;
		string name("output_files/model_");
		string suf(".pdb");
		ss << i;
		ss >> str;
		Enhance_population[vec_dscore_i[i].second].dump_pdb(name+str+suf);
	}
}

void 
ClassicAbinitio::Enhance(){
	core::pose::Pose pose_E;
	dscore_map.resize( lenth_of_sequence );
	for(int i = 0; i < lenth_of_sequence; i++){
		for(int j = 0; j < lenth_of_sequence; j++){
			dscore_map[i].push_back(0);
		}
	}
	for(int i = 0; i < Enhance_population.size(); i++){
		for(int position_cycle_in = 1; position_cycle_in <= 1; position_cycle_in++){
			vector<Size> shib_po;
			shib_po.push_back(0);
			for(int position_cycle = 0; position_cycle < lenth_of_sequence; position_cycle++){
				bool cheng(false);
				
				lock_position(shib_po);
				
				double last_max_vec_mean = max_vec_mean;
				FragAssem_ = brute_move_large_;
				
				for(int fragment_cycle = 0; fragment_cycle < 30; fragment_cycle++){
					if(max_vec_mean >= last_max_vec_mean ){
						if(max_vec_mean_index+1>0 && max_vec_mean_index+1<lenth_of_sequence-8){
							core::pose::Pose pose_f;
							pose_f = Enhance_population[i];
							pose_E = Enhance_population[i];
							
							bml->define_start_window(max_vec_mean_index+1);
							FragAssem_->apply( pose_f );
							
							Real temp1 = pose_f.phi(max_vec_mean_index+1);
							Real temp2 = pose_f.psi(max_vec_mean_index+1);
							Real temp3 = pose_f.omega(max_vec_mean_index+1);
							pose_E.set_phi(max_vec_mean_index+1,temp1);
							pose_E.set_psi(max_vec_mean_index+1,temp2);
							pose_E.set_omega(max_vec_mean_index+1,temp3);
							
							Enhance_Distance_score(pose_E);
							max_vec_mean = again_calc();
						}
					}
					else{
						Enhance_population[i] = pose_E;
						cheng = true;
						
						break;
					}
				}
				
				if(!cheng){
					shib_po.push_back(max_vec_mean_index);
				}
			}
		}
	  
	  
		for(int position_cycle_in = 1; position_cycle_in <= 1; position_cycle_in++){
			vector<Size> shib_po;
			shib_po.push_back(0);
			for(int position_cycle = 0; position_cycle < lenth_of_sequence; position_cycle++){
				bool cheng(false);
			
				lock_position(shib_po);
				
				double last_max_vec_mean = max_vec_mean;
				FragAssem_ = brute_move_small_;
				
				for(int fragment_cycle = 0; fragment_cycle < 30; fragment_cycle++){
					if(max_vec_mean >= last_max_vec_mean ){
						core::pose::Pose pose_f;
						pose_f = Enhance_population[i];
						pose_E = Enhance_population[i];
						
						bms->define_start_window(max_vec_mean_index+1);
						FragAssem_->apply( pose_f );
						
						Real temp1 = pose_f.phi(max_vec_mean_index+1);
						Real temp2 = pose_f.psi(max_vec_mean_index+1);
						Real temp3 = pose_f.omega(max_vec_mean_index+1);
						pose_E.set_phi(max_vec_mean_index+1,temp1);
						pose_E.set_psi(max_vec_mean_index+1,temp2);
						pose_E.set_omega(max_vec_mean_index+1,temp3);
						
						Enhance_Distance_score(pose_E);
						max_vec_mean = again_calc();
					}
					else{
						Enhance_population[i] = pose_E;
						cheng = true;
						
						break;
					}
				}
				
				if(!cheng){
					shib_po.push_back(max_vec_mean_index);
				}
			}
		}
		
		for(int position_cycle_in = 1; position_cycle_in <= 1; position_cycle_in++){
			vector<Size> shib_po;
			shib_po.push_back(0);
			for(int position_cycle = 0; position_cycle < lenth_of_sequence; position_cycle++){
				bool cheng(false);
				
				lock_position(shib_po);
				
				double last_max_vec_mean = max_vec_mean;
				FragAssem_ = brute_move_small_;
				
				for(int fragment_cycle = 0; fragment_cycle < 10; fragment_cycle++){
					if(max_vec_mean >= last_max_vec_mean ){
						core::pose::Pose pose_fa;
						core::pose::Pose pose_fb;
						pose_fa = Enhance_population[i];
						pose_fb = Enhance_population[i];
						pose_E = Enhance_population[i];
						
						bms->define_start_window(max_vec_mean_index+1);
						FragAssem_->apply( pose_fa );
						FragAssem_->apply( pose_fb );
						
						Real temp1 = pose_E.phi(max_vec_mean_index+1) + 0.5*( pose_fb.phi(max_vec_mean_index+1)-pose_fa.phi(max_vec_mean_index+1) );
						Real temp2 = pose_E.psi(max_vec_mean_index+1) + 0.5*( pose_fb.psi(max_vec_mean_index+1)-pose_fa.psi(max_vec_mean_index+1) );
						Real temp3 = pose_E.omega(max_vec_mean_index+1) + 0.5*( pose_fb.omega(max_vec_mean_index+1)-pose_fa.omega(max_vec_mean_index+1) );
						pose_E.set_phi(max_vec_mean_index+1,temp1);
						pose_E.set_psi(max_vec_mean_index+1,temp2);
						pose_E.set_omega(max_vec_mean_index+1,temp3);
						
						Enhance_Distance_score(pose_E);
						max_vec_mean = again_calc();
					}
					else{
						Enhance_population[i] = pose_E;
						cheng = true;
						
						break;
					}
				}
				
				if(!cheng){
					shib_po.push_back(max_vec_mean_index);
				}
			}
		}
		
		for(int position_cycle_in = 1; position_cycle_in <= 1; position_cycle_in++){
			vector<Size> shib_po;
			shib_po.push_back(0);
			for(int position_cycle = 0; position_cycle < lenth_of_sequence; position_cycle++){
				bool cheng(false);
				
				lock_position(shib_po);
				
				double last_max_vec_mean = max_vec_mean;
				FragAssem_ = brute_move_small_;
			
				for(int fragment_cycle = 0; fragment_cycle < 10; fragment_cycle++){
					if(max_vec_mean >= last_max_vec_mean ){
						core::pose::Pose pose_fa;
						core::pose::Pose pose_fb;
						pose_fa = Enhance_population[i];
						pose_fb = Enhance_population[i];
						pose_E = Enhance_population[i];
						
						bms->define_start_window(max_vec_mean_index+1);
						FragAssem_->apply( pose_fa );
						FragAssem_->apply( pose_fb );
						
						Real temp1 = pose_E.phi(max_vec_mean_index+1) + 0.5*( pose_E.phi(max_vec_mean_index+1)-pose_fa.phi(max_vec_mean_index+1) );
						Real temp2 = pose_E.psi(max_vec_mean_index+1) + 0.5*( pose_E.psi(max_vec_mean_index+1)-pose_fa.psi(max_vec_mean_index+1) );
						Real temp3 = pose_E.omega(max_vec_mean_index+1) + 0.5*( pose_E.omega(max_vec_mean_index+1)-pose_fa.omega(max_vec_mean_index+1) );
						pose_E.set_phi(max_vec_mean_index+1,temp1);
						pose_E.set_psi(max_vec_mean_index+1,temp2);
						pose_E.set_omega(max_vec_mean_index+1,temp3);
						
						Enhance_Distance_score(pose_E);
						max_vec_mean = again_calc();
					}
					else{
						Enhance_population[i] = pose_E;
						cheng = true;
						
						break;
					}
				}
				
				if(!cheng){
					shib_po.push_back(max_vec_mean_index);
				}
			}
		}
	}
}

void 
ClassicAbinitio::lock_position(vector<Size> shib){
	vector<double> vec_mean;
	for(int i = 0; i < lenth_of_sequence; i++){
		if( i == 0 || i>=lenth_of_sequence-3 ){
			vec_mean.push_back(0);
		}
		else{
			int l =0;
			double local_score = 0;
			for(int m = i; m < lenth_of_sequence; m++){
				for(int n = 0; n < i; n++){
					if( dscore_map[m][n] !=0 ){
						local_score += dscore_map[m][n];
						l++;
					}
				}
			}
			double mean_local_score = 0;
			if( l!=0 ){
				mean_local_score = local_score/l;
			}
			vec_mean.push_back(mean_local_score);
		}
	}
	for(Size i = 0; i < shib.size(); i++){
		vec_mean[shib[i]] = 0;
	}
	
	max_vec_mean_index = max_element(vec_mean.begin(), vec_mean.end())-vec_mean.begin();
	max_vec_mean = *max_element(vec_mean.begin(), vec_mean.end());
}

double 
ClassicAbinitio::again_calc(){
	int l = 0;
	double local_score = 0;
	for(int m = max_vec_mean_index; m < lenth_of_sequence; m++){
		for(Size n = 0; n < max_vec_mean_index; n++){
			if( dscore_map[m][n] !=0 ){
				local_score += dscore_map[m][n];
				l++;
			}
		}
	}
	
	double mean_local_score = 0;
	if( l!=0 ){
		mean_local_score = local_score/l;
	}
	
	return mean_local_score;
}

double 
ClassicAbinitio::Enhance_Distance_score(core::pose::Pose& pose){
	double Dis_score = 0;
	double confidence_all = 0;
	double d0 = 0;
	double di = 0;
	for(Size k = 0; k < vec_distance.size(); k++){
		double C_distance = 0;
		double ijscore = 0;
		if ( pose.residue(vec_distance[k].re1).name3() == "GLY" || pose.residue(vec_distance[k].re2).name3() == "GLY" ){
			numeric::xyzVector<Real> CA1 = pose.residue(vec_distance[k].re1).xyz("CA");
			numeric::xyzVector<Real> CA2 = pose.residue(vec_distance[k].re2).xyz("CA");
			C_distance = CA1.distance(CA2);
		}else{
			numeric::xyzVector<Real> CB1 = pose.residue(vec_distance[k].re1).xyz("CB");
			numeric::xyzVector<Real> CB2 = pose.residue(vec_distance[k].re2).xyz("CB");
			C_distance = CB1.distance(CB2);
		}
		if( fabs(vec_distance[k].re1 - vec_distance[k].re2) >= 3 ){
			d0 = log(fabs(vec_distance[k].re1 - vec_distance[k].re2));
			di = fabs(C_distance - vec_distance[k].dist);
			ijscore = vec_distance[k].confidence*(di/d0);
			Dis_score += ijscore;
			dscore_map[vec_distance[k].re1 - 1][vec_distance[k].re2 - 1] = ijscore;
		}
	}
	return Dis_score;
}

void 
ClassicAbinitio::Read_parameters(){
	ifstream Rp("./parameter");
	string line;
	while (getline(Rp,line)){
		if ( line[0] != '#' && line[0] != ' ' ){
			istringstream line_data( line );
			string parameterName;
			int parameterValue;
			line_data >> parameterName >> parameterValue;
			parametersMap.insert(make_pair(parameterName, parameterValue));
		}
	}
	Rp.close();
}

void 
ClassicAbinitio::get_parameters(map<string, int>& parametersMap){
	if (parametersMap.find("N") != parametersMap.end())
		N = parametersMap["N"];
	else{
		cout << "====================================parameter 'N' is not found in parameters!!!" << endl;
		exit(0);
	}
	
	if (parametersMap.find("Crowding_radius_exp") != parametersMap.end())
	    Crowding_radius_exp = parametersMap["Crowding_radius_exp"];
	else{
	    cout << "====================================parameter 'Crowding_radius_exp' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	if (parametersMap.find("Crowding_radius_mer") != parametersMap.end())
	    Crowding_radius_mer = parametersMap["Crowding_radius_mer"];
	else{
	    cout << "====================================parameter 'Crowding_radius_mer' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	if (parametersMap.find("G_explore") != parametersMap.end())
	    G_explore = parametersMap["G_explore"];
	else{
	    cout << "====================================parameter 'G_explore' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	if (parametersMap.find("G_maintain") != parametersMap.end())
	    G_maintain = parametersMap["G_maintain"];
	else{
	    cout << "====================================parameter 'G_maintain' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	if (parametersMap.find("G_enhance_c") != parametersMap.end())
	    G_enhance_c = parametersMap["G_enhance_c"];
	else{
	    cout << "====================================parameter 'G_enhance_c' is not found in parameters!!!" << endl;
	    exit(0);
	}
}


vector< core::pose::Pose >
ClassicAbinitio::Initial_population(core::pose::Pose &pose){
	pose::Pose ready_pose = pose;
	for(int i = 0; i < N; i++){
		RMA_stage1( pose );
		pose = ready_pose;
	}
	
	vector< core::pose::Pose > storage_pose;
	RMA_stage2( pose );
	for (int i = 0; i < N; i++){
		storage_pose.push_back(All_pose[i]);
	}
	
	return storage_pose;
}


void ClassicAbinitio::RMA_stage1( core::pose::Pose &pose ) {
	using namespace moves;
	using namespace scoring;
	using namespace scoring::constraints;
	bool success( true );
	if ( !bSkipStage1_ ) {
		PROF_START( basic::STAGE1 );
		clock_t starttime = clock();

		if ( !prepare_stage1( pose ) ) {
			set_last_move_status( moves::FAIL_RETRY );
			return;
		}
		// part 1 ----------------------------------------
		tr.Info <<  "\n===================================================================\n";
		tr.Info <<  "   Stage 1                                                         \n";
		tr.Info <<  "   Folding with score0 for max of " << stage1_cycles() << std::endl;

		if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
		if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
			output_debug_structure( pose, "stage0" );
		}
		if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(), "stage_1", false /* fullatom*/, true /*fold tree */ ) ) {
			ConstraintSetOP orig_constraints(nullptr);
			orig_constraints = pose.constraint_set()->clone();
			success = do_stage1_cycles( pose );

			if ( tr.Info.visible() ) current_scorefxn().show( tr, pose );
			recover_low( pose, STAGE_1 );
			mc().show_counters();
			total_trials_+=mc().total_trials();
			mc().reset_counters();

			pose.constraint_set( orig_constraints ); // restore constraints - this is critical for checkpointing to work
			get_checkpoints().checkpoint( pose, get_current_tag(), "stage_1", true /*fold tree */ );
		} //recover checkpoint
		get_checkpoints().debug( get_current_tag(), "stage_1", current_scorefxn()( pose ) );

		clock_t endtime = clock();
		PROF_STOP( basic::STAGE1 );
		if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
		if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
			tr.Info << "Timeperstep: " << (double(endtime) - starttime )/(CLOCKS_PER_SEC ) << std::endl;
			output_debug_structure( pose, "stage1" );
		}
	} //skipStage1
	if ( !success ) {
		set_last_move_status( moves::FAIL_RETRY );
		return;
	}
	All_pose.push_back(pose);
}

void ClassicAbinitio::RMA_stage2(core::pose::Pose &pose ) {
	using namespace moves;
	using namespace scoring;
	using namespace scoring::constraints;
	bool success( true );
	if ( !bSkipStage2_ ) {
		for(int i=0; i<N; i++){
			//
			// part 2 ----------------------------------------
			tr.Info <<  "\n===================================================================\n";
			tr.Info <<  "   Stage 2" << "   N = " << i << endl;
			tr.Info <<  "   Folding with score1 for " << stage2_cycles() << std::endl;

			PROF_START( basic::STAGE2 );
			clock_t starttime = clock();
			
			pose = All_pose[i];

			if ( close_chbrk_ ) {
				Real const setting( 0.25 );
				set_score_weight( scoring::linear_chainbreak, setting, STAGE_2 );
				tr.Info <<  " Chain_break score assigned " << std::endl;
			}


			if ( !prepare_stage2( pose ) )  {
				set_last_move_status( moves::FAIL_RETRY );
				return;
			}

			if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(), "stage_2", false /* fullatom */, true /*fold tree */ ) ) {
				ConstraintSetOP orig_constraints(nullptr);
				orig_constraints = pose.constraint_set()->clone();

				success = do_stage2_cycles( pose );
				recover_low( pose, STAGE_2 );                    //default OFF: seems to be a bad choice after score0

				if  ( tr.visible() ) current_scorefxn().show( tr, pose );
				mc().show_counters();
				total_trials_+=mc().total_trials();
				mc().reset_counters();

				pose.constraint_set( orig_constraints ); // restore constraints - this is critical for checkpointing to work
				get_checkpoints().checkpoint( pose, get_current_tag(), "stage_2", true /*fold tree */ );
			}
			get_checkpoints().debug( get_current_tag(), "stage_2", current_scorefxn()( pose ) );

			clock_t endtime = clock();
			PROF_STOP( basic::STAGE2 );
			if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
			if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
				output_debug_structure( pose, "stage2" );
				tr << "Timeperstep: " << (double(endtime) - starttime )/(CLOCKS_PER_SEC ) << std::endl;
			}
			
			All_pose[i]=pose;
		}
	} //bSkipStage2

	if ( !success ) {
		set_last_move_status( moves::FAIL_RETRY );
		return;
	}
}

void ClassicAbinitio::RMA_stage3( core::pose::Pose &pose ) {
	using namespace moves;
	using namespace scoring;
	using namespace scoring::constraints;
	bool success( true );
	if ( !bSkipStage3_ ) {
		// moved checkpointing into do_stage3_cycles because of structure store
		for(int i=0;i<N;i++){
			// part 3 ----------------------------------------
			tr.Info <<  "\n===================================================================\n";
			tr.Info <<  "   Stage 3" << "   N = " << i << endl;
			tr.Info <<  "   Folding with score2 and score5 for " << stage3_cycles() <<std::endl;

			PROF_START( basic::STAGE3 );
			clock_t starttime = clock();
			
			pose = All_pose[i];

			if ( !prepare_stage3( pose ) ) {
				set_last_move_status( moves::FAIL_RETRY );
				return;
			}
			// this is not the final score-function.. only known after prepare_loop_in_stage3
			// because this is confusing rather not show.if ( tr.Info.visible() ) current_scorefxn().show( tr, pose );

			success = do_stage3_cycles( pose );
			recover_low( pose, STAGE_3b );

			if ( tr.Info.visible() ) current_scorefxn().show( tr, pose );
			mc().show_counters();
			total_trials_+=mc().total_trials();
			mc().reset_counters();

			clock_t endtime = clock();
			PROF_STOP( basic::STAGE3);
			if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
			if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
				output_debug_structure( pose, "stage3" );
				tr << "Timeperstep: " << (double(endtime) - starttime )/( CLOCKS_PER_SEC) << std::endl;
			}

			//  pose.dump_pdb("stage3.pdb");
			
			All_pose[i]=pose;
		}

	}

	if ( !success ) {
		set_last_move_status( moves::FAIL_RETRY );
		return;
	}
}

void ClassicAbinitio::RMA_stage4( core::pose::Pose &pose ) {
	using namespace moves;
	using namespace scoring;
	using namespace scoring::constraints;
	if ( !bSkipStage4_ ) {		
		for(int i=0;i<N;i++){

			// part 4 ------------------------------------------
			tr.Info <<  "\n===================================================================\n";
			tr.Info <<  "   Stage 4" << "   N = " << i << endl;
			tr.Info <<  "   Folding with score3 for " << stage4_cycles() <<std::endl;

			PROF_START( basic::STAGE4 );
			clock_t starttime = clock();

			if ( !prepare_stage4( pose ) ) {
				set_last_move_status( moves::FAIL_RETRY );
				return;
			}
			
			pose = All_pose[i];

			//score-fxn may be changed in do_stage4_cycles...
			// confusing if shown here already... if ( tr.Info.visible() ) current_scorefxn().show( tr, pose);
			bool success = do_stage4_cycles( pose );
			recover_low( pose, STAGE_4  );

			if ( tr.Info.visible() ) current_scorefxn().show( tr, pose);
			mc().show_counters();
			total_trials_+=mc().total_trials();
			mc().reset_counters();

			clock_t endtime = clock();
			PROF_STOP( basic::STAGE4 );
			if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
			if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
				output_debug_structure( pose, "stage4" );
				tr << "Timeperstep: " << (double(endtime) - starttime )/( CLOCKS_PER_SEC ) << std::endl;
			}

			tr.Info <<  "\n===================================================================\n";
			tr.Info <<  "   Finished Abinitio                                                 \n";
			tr.Info <<  std::endl;
			//  pose.dump_pdb("stage4.pdb");
			
			All_pose[i]=pose;
		}
	}
}

void ClassicAbinitio::RMA_stage5( core::pose::Pose &pose ) {
	using namespace moves;
	using namespace scoring;
	using namespace scoring::constraints;
	bool success( true );
	if ( !bSkipStage5_ ) {

		// part 5 ------------------------------------------
		tr.Info <<  "\n===================================================================\n";
		tr.Info <<  "   Stage 5                                                         \n";
		tr.Info <<  "   Folding with score3 for " << stage5_cycles() <<std::endl;

		PROF_START( basic::STAGE5 );
		clock_t starttime = clock();

		if ( !prepare_stage5( pose ) ) {
			set_last_move_status( moves::FAIL_RETRY );
			return;
		}

		success = do_stage5_cycles( pose );
		recover_low( pose, STAGE_5 );

		if ( tr.Info.visible() ) current_scorefxn().show( tr, pose);
		//  current_scorefxn().show(tr, pose);
		mc().show_counters();
		total_trials_+=mc().total_trials();
		mc().reset_counters();

		clock_t endtime = clock();
		PROF_STOP( basic::STAGE5 );
		if ( option[ basic::options::OptionKeys::run::profile ] ) prof_show();
		if ( option[ basic::options::OptionKeys::abinitio::debug ]() ) {
			output_debug_structure( pose, "stage5" );
			tr << "Timeperstep: " << (double(endtime) - starttime )/( CLOCKS_PER_SEC ) << std::endl;
		}

		tr.Info <<  "\n===================================================================\n";
		tr.Info <<  "   Now really finished Abinitio                                                 \n";
		tr.Info <<  std::endl;
		//  pose.dump_pdb("stage5.pdb");
	}

	get_checkpoints().flush_checkpoints();
	if ( !success ) set_last_move_status( moves::FAIL_RETRY );
	//basic::prof_show();

	return;
}// ClassicAbinitio::apply( pose::Pose & pose )


std::string
ClassicAbinitio::get_name() const {
	return "ClassicAbinitio";
}

//@brief return FramgentMover for smooth_small fragment insertions (i.e., stage4 moves)
simple_moves::FragmentMoverOP
ClassicAbinitio::smooth_move_small() {
	return smooth_move_small_;
}

//@brief return FragmentMover for small fragment insertions ( i.e., stage3/4 moves )
simple_moves::FragmentMoverOP
ClassicAbinitio::brute_move_small() {
	return brute_move_small_;
}

//@brief return FragmentMover for large fragment insertions (i.e., stage1/2 moves )
simple_moves::FragmentMoverOP
ClassicAbinitio::brute_move_large() {
	return brute_move_large_;
}

//@brief change the movemap ( is propagated to mover-objects )
//@detail overload if your extension stores additional moves as member variables
void
ClassicAbinitio::set_movemap( core::kinematics::MoveMapCOP mm )
{
	movemap_ = mm;
	if ( smooth_move_small_ ) smooth_move_small_->set_movemap( mm );
	if ( brute_move_small_  ) brute_move_small_ ->set_movemap( mm );
	if ( brute_move_large_  ) brute_move_large_ ->set_movemap( mm );
}

//@brief set new instances of FragmentMovers
void
ClassicAbinitio::set_moves(
	simple_moves::FragmentMoverOP brute_move_small,
	simple_moves::FragmentMoverOP brute_move_large,
	simple_moves::FragmentMoverOP smooth_move_small
)
{
	smooth_move_small_ = smooth_move_small;
	brute_move_small_  = brute_move_small;
	brute_move_large_  = brute_move_large;
	update_moves();
}

//@brief returns current movemap
core::kinematics::MoveMapCOP
ClassicAbinitio::movemap() {
	return movemap_;
}

//@detail read cmd_line options and set default versions for many protocol members: trials/moves, score-functions, Monte-Carlo
void ClassicAbinitio::set_defaults( pose::Pose const& pose ) {
	temperature_ = 2.0;
	bSkipStage1_ = false;
	bSkipStage2_ = false;
	bSkipStage3_ = false;
	bSkipStage4_ = false;
	bSkipStage5_ = true; //vats is turned off by default
	set_default_scores();
	set_default_options();
	set_default_mc( pose, *score_stage1_ );
	update_moves();
}

//@detail called to notify about changes in Movers: new movemap or Moverclass
void ClassicAbinitio::update_moves() {
	/* set apply_large_frags_ and
	short_insert_region_
	*/
	/* what about move-map ? It can be set manually for all Fragment_Moves .. */
	// set_move_map();
	set_trials();
}

//@detail create instances of TrialMover for our FragmentMover objects
void ClassicAbinitio::set_trials() {
	// setup loop1
	runtime_assert( brute_move_large_ != nullptr );
	trial_large_ = moves::TrialMoverOP( new moves::TrialMover( brute_move_large_, mc_ ) );
	//trial_large_->set_keep_stats( true );
	trial_large_->keep_stats_type( moves::accept_reject );

	runtime_assert( brute_move_small_ != nullptr );
	trial_small_ = moves::TrialMoverOP( new moves::TrialMover( brute_move_small_, mc_ ) );
	//trial_small_->set_keep_stats( true );
	trial_small_->keep_stats_type( moves::accept_reject );

	runtime_assert( smooth_move_small_ != nullptr );
	smooth_trial_small_ = moves::TrialMoverOP( new moves::TrialMover( smooth_move_small_, mc_ ) );
	//smooth_trial_small_->set_keep_stats( true );
	smooth_trial_small_->keep_stats_type( moves::accept_reject );

	//build trial_pack mover
	moves::SequenceMoverOP combo_small( new moves::SequenceMover() );
	combo_small->add_mover(brute_move_small_);
	combo_small->add_mover(pack_rotamers_);
	trial_small_pack_ = moves::TrialMoverOP( new moves::TrialMover(combo_small, mc_) );
	moves::SequenceMoverOP combo_smooth( new moves::SequenceMover() );
	combo_smooth->add_mover(smooth_move_small_);
	combo_smooth->add_mover(pack_rotamers_);
	smooth_trial_small_pack_ = moves::TrialMoverOP( new moves::TrialMover(combo_smooth, mc_) );
}

//@detail sets Monto-Carlo object to default
void ClassicAbinitio::set_default_mc(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn
) {
	set_mc( moves::MonteCarloOP( new moves::MonteCarlo( pose, scorefxn, temperature_ ) ) );
}

//@detail sets Monto-Carlo object
void ClassicAbinitio::set_mc( moves::MonteCarloOP mc_in ) {
	mc_ = mc_in;
	if ( trial_large_ ) trial_large_->set_mc( mc_ );
	if ( trial_small_ ) trial_small_->set_mc( mc_ );
	if ( smooth_trial_small_ ) smooth_trial_small_->set_mc( mc_ );
}

//@detail override cmd-line setting for "increase_cycling"
void ClassicAbinitio::set_cycles( Real increase_cycles ) {
	stage1_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage2_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage3_cycles_ = static_cast< int > (2000 * increase_cycles);
	stage4_cycles_ = static_cast< int > (4000 * increase_cycles);
	stage5_cycles_ = static_cast< int > (50000* increase_cycles);//vats

	using namespace basic::options;
	if ( option[ OptionKeys::abinitio::only_stage1 ]() ) {
		stage2_cycles_ = 0;
		stage3_cycles_ = 0;
		stage4_cycles_ = 0;
		bSkipStage2_ = bSkipStage3_ = /*bSkipStage3_ =*/ true;  // Was bSkipStage4_ meant? ~Labonte
	}
}

void ClassicAbinitio::set_default_scores() {
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	tr.Debug << "creating standard scoring functions" << std::endl;

	if ( option[ OptionKeys::abinitio::stage1_patch ].user() ) {
		score_stage1_  = ScoreFunctionFactory::create_score_function( "score0", option[ OptionKeys::abinitio::stage1_patch ]() );
	} else {
		score_stage1_  = ScoreFunctionFactory::create_score_function( "score0" );
	}

	if ( option[ OptionKeys::abinitio::stage2_patch ].user() ) {
		score_stage2_  = ScoreFunctionFactory::create_score_function( "score1", option[ OptionKeys::abinitio::stage2_patch ]() );
	} else {
		score_stage2_  = ScoreFunctionFactory::create_score_function( "score1" );
	}

	if ( option[ OptionKeys::abinitio::stage3a_patch ].user() ) {
		score_stage3a_ = ScoreFunctionFactory::create_score_function( "score2", option[ OptionKeys::abinitio::stage3a_patch ]() );
	} else {
		score_stage3a_ = ScoreFunctionFactory::create_score_function( "score2" );
	}

	if ( option[ OptionKeys::abinitio::stage3b_patch ].user() ) {
		score_stage3b_ = ScoreFunctionFactory::create_score_function( "score5", option[ OptionKeys::abinitio::stage3b_patch ]() );
	} else {
		score_stage3b_ = ScoreFunctionFactory::create_score_function( "score5" );
	}

	if ( option[ OptionKeys::abinitio::stage4_patch ].user() ) {
		score_stage4_  = ScoreFunctionFactory::create_score_function( "score3", option[ OptionKeys::abinitio::stage4_patch ]() );
	} else {
		score_stage4_  = ScoreFunctionFactory::create_score_function( "score3" );
	}

	//loading the cenrot score
	score_stage4rot_ = ScoreFunctionFactory::create_score_function( "score4_cenrot_relax" );
	//score_stage4rot_->set_weight(core::scoring::cen_rot_dun, 0.0);
	score_stage4rot_sc_ = ScoreFunctionFactory::create_score_function( "score4_cenrot_repack" );
	//score_stage4rot_sc_->set_weight(core::scoring::cen_rot_dun, 1.0);

	if ( option[ OptionKeys::abinitio::stage5_patch ].user() ) { //vats
		score_stage5_  = ScoreFunctionFactory::create_score_function( "score3", option[ OptionKeys::abinitio::stage5_patch ]() );
	} else {
		score_stage5_  = ScoreFunctionFactory::create_score_function( "score3" );
	}


	if ( option[ OptionKeys::abinitio::override_vdw_all_stages ] ) {
		set_score_weight( scoring::vdw, option[ OptionKeys::abinitio::vdw_weight_stage1 ], ALL_STAGES );
	}
}


/// @brief sets a score weight for all stages of abinitio
void ClassicAbinitio::set_score_weight( scoring::ScoreType type, Real setting, StageID stage ) {
	tr.Debug << "set score weights for ";
	if ( stage == ALL_STAGES ) tr.Debug << "all stages ";
	else tr.Debug << "stage " << (stage <= STAGE_3a ? stage : ( stage-1 ) ) << ( stage == STAGE_3b ? "b " : " " );
	tr.Debug << scoring::name_from_score_type(type) << " " << setting << std::endl;
	if ( score_stage1_  && ( stage == STAGE_1  || stage == ALL_STAGES ) ) score_stage1_ ->set_weight(type, setting);
	if ( score_stage2_  && ( stage == STAGE_2  || stage == ALL_STAGES ) ) score_stage2_ ->set_weight(type, setting);
	if ( score_stage3a_ && ( stage == STAGE_3a || stage == ALL_STAGES ) ) score_stage3a_->set_weight(type, setting);
	if ( score_stage3b_ && ( stage == STAGE_3b || stage == ALL_STAGES ) ) score_stage3b_->set_weight(type, setting);
	if ( score_stage4_  && ( stage == STAGE_4  || stage == ALL_STAGES ) ) score_stage4_ ->set_weight(type, setting);
	if ( score_stage4rot_  && ( stage == STAGE_4  || stage == ALL_STAGES ) ) score_stage4rot_ ->set_weight(type, setting);
	if ( score_stage5_  && ( stage == STAGE_5  || stage == ALL_STAGES ) ) score_stage5_ ->set_weight(type, setting);//vats
}

//@brief currently used score function ( depends on stage )
scoring::ScoreFunction const& ClassicAbinitio::current_scorefxn() const {
	return mc().score_function();
}

//@brief set current scorefunction
void ClassicAbinitio::current_scorefxn( scoring::ScoreFunction const& scorefxn ) {
	mc().score_function( scorefxn );
}

//@brief set individual weight of current scorefunction --- does not change the predefined scores: score_stageX_
void ClassicAbinitio::set_current_weight( core::scoring::ScoreType type, core::Real setting ) {
	scoring::ScoreFunctionOP scorefxn ( mc().score_function().clone() );
	scorefxn->set_weight( type, setting );
	mc().score_function( *scorefxn ); //trigger rescore
}

void ClassicAbinitio::set_default_options() {
	bSkipStage1_ = bSkipStage2_ = bSkipStage3_ = bSkipStage4_ = false;
	bSkipStage5_ = true; //vats turned off by default
	using namespace basic::options;
	just_smooth_cycles_ = option[ OptionKeys::abinitio::smooth_cycles_only ]; // defaults to false
	bQuickTest_ = basic::options::option[ basic::options::OptionKeys::run::test_cycles ]();

	if ( bQuickTest() ) {
		set_cycles( 0.001 );
	} else {
		set_cycles( option[ OptionKeys::abinitio::increase_cycles ] ); // defaults to factor of 1.0
	}

	if ( just_smooth_cycles_ ) {
		bSkipStage1_ = bSkipStage2_ = bSkipStage3_ = bSkipStage5_ = true;
	}
	if ( option[ OptionKeys::abinitio::only_stage1 ] ) {
		bSkipStage2_ = bSkipStage3_ = bSkipStage4_ = bSkipStage5_= true;
	}

	if ( option[ OptionKeys::abinitio::include_stage5 ] ) {
		bSkipStage5_ = false;
	}

	apply_large_frags_   = true;  // apply large frags in phase 2!

	// in rosetta++ switched on in fold_abinitio if contig_size < 30 in pose_abinitio never
	short_insert_region_ = false;  // apply small fragments in phase 2!

	if ( option[ OptionKeys::abinitio::recover_low_in_stages ].user() ) {
		for ( int it : option[ OptionKeys::abinitio::recover_low_in_stages ]() ) {
			if ( it == 1 ) recover_low_stages_.push_back( STAGE_1 );
			else if ( it == 2 ) recover_low_stages_.push_back( STAGE_2 );
			else if ( it == 3 ) {
				recover_low_stages_.push_back( STAGE_3a );
				recover_low_stages_.push_back( STAGE_3b );
			} else if ( it == 4 ) recover_low_stages_.push_back( STAGE_4 );
		}
	} else {
		recover_low_stages_.clear();
		recover_low_stages_.push_back( STAGE_1 );
		recover_low_stages_.push_back( STAGE_2 );
		recover_low_stages_.push_back( STAGE_3a );
		recover_low_stages_.push_back( STAGE_3b );
		recover_low_stages_.push_back( STAGE_4 );
		recover_low_stages_.push_back( STAGE_5 );
	}

	close_chbrk_ = option[ OptionKeys::abinitio::close_chbrk ];

}


/// @brief (helper) functor class which keeps track of old pose for the
/// convergence check in stage3 cycles
/// @detail
/// calls of operator ( pose ) compare the
class hConvergenceCheck;
using hConvergenceCheckOP = utility::pointer::shared_ptr<hConvergenceCheck>;

class hConvergenceCheck : public moves::PoseCondition {
public:
	hConvergenceCheck()= default;
	void reset() { ct_ = 0; bInit_ = false; }
	void set_trials( moves::TrialMoverOP trin ) {
		trials_ = trin;
		runtime_assert( trials_->keep_stats_type() < moves::no_stats );
		last_move_ = 0;
	}
	bool operator() ( const core::pose::Pose & pose ) override;
private:
	pose::Pose very_old_pose_;
	bool bInit_{ false };
	Size ct_ = 0;
	moves::TrialMoverOP trials_;
	Size last_move_;
};

// keep going --> return true
bool hConvergenceCheck::operator() ( const core::pose::Pose & pose ) {
	if ( !bInit_ ) {
		bInit_ = true;
		very_old_pose_ = pose;
		return true;
	}
	runtime_assert( trials_ != nullptr );
	tr.Trace << "TrialCounter in hConvergenceCheck: " << trials_->num_accepts() << std::endl;
	if ( numeric::mod(trials_->num_accepts(),100) != 0 ) return true;
	if ( (Size) trials_->num_accepts() <= last_move_ ) return true;
	last_move_ = trials_->num_accepts();
	// change this later to this: (after we compared with rosetta++ and are happy)
	// if ( numeric::mod(++ct_, 1000) != 0 ) return false; //assumes an approx acceptance rate of 0.1

	// still here? do the check:

	core::Real converge_rms = core::scoring::CA_rmsd( very_old_pose_, pose );
	very_old_pose_ = pose;
	if ( converge_rms >= 3.0 ) {
		return true;
	}
	// if we get here thing is converged stop the While-Loop
	tr.Info << " stop cycles in stage3 due to convergence " << std::endl;
	return false;
}


bool ClassicAbinitio::do_stage1_cycles( pose::Pose &pose ) {
	AllResiduesChanged done( pose, brute_move_large()->insert_map(), *movemap() );
	moves::MoverOP trial( stage1_mover( pose, trial_large() ) );

	// FragmentMoverOP frag_mover = brute_move_large_;
	// fragment::FragmentIO().write("stage1_frags_classic.dat",*frag_mover->fragments());

	Size j;
	for ( j = 1; j <= stage1_cycles(); ++j ) {
		trial->apply( pose ); // apply a large fragment insertion, accept with MC boltzmann probability
		if ( done(pose) ) {
			tr.Info << "Replaced extended chain after " << j << " cycles." << std::endl;
			mc().reset( pose ); // make sure that we keep the final structure
			return true;
		}
	}
	tr.Warning << "extended chain may still remain after " << stage1_cycles() << " cycles!" << std::endl;
	done.show_unmoved( pose, tr.Warning );
	mc().reset( pose ); // make sure that we keep the final structure
	return true;
}

bool ClassicAbinitio::do_stage2_cycles( pose::Pose &pose ) {

	//setup cycle
	moves::SequenceMoverOP cycle( new moves::SequenceMover() );
	if ( apply_large_frags_   ) cycle->add_mover( trial_large_->mover() );
	if ( short_insert_region_ ) cycle->add_mover( trial_small_->mover() );

	Size nr_cycles = stage2_cycles() / ( short_insert_region_ ? 2 : 1 );
	moves::TrialMoverOP trials( new moves::TrialMover( cycle, mc_ptr() ) );
	moves::RepeatMover( stage2_mover( pose, trials ), nr_cycles ).apply(pose);

	//is there a better way to find out how many steps ? for instance how many calls to scoring?
	return true; // as best guess
}

/*! @detail stage3 cycles:
nloop1 : outer iterations
nloop2 : inner iterations
stage3_cycle : trials per inner iteration
every inner iteration we switch between score_stage3a ( default: score2 ) and score_stage3b ( default: score 5 )

prepare_loop_in_stage3() is called before the stage3_cycles() of trials are started.

first outer loop-iteration is done with TrialMover trial_large()
all following iterations with trial_small()

start each iteration with the lowest_score_pose. ( mc->recover_low() -- called in prepare_loop_in_stage3() )

*/
bool ClassicAbinitio::do_stage3_cycles( pose::Pose &pose ) {
	using namespace ObjexxFCL;

	// interlaced score2 / score 5 loops
	// nloops1 and nloops2 could become member-variables and thus changeable from the outside
	int nloop1 = 1;
	int nloop2 = 10; //careful: if you change these the number of structures in the structure store changes.. problem with checkpointing
	// individual checkpoints for each stage3 iteration would be a remedy. ...

	if ( short_insert_region_ ) {
		nloop1 = 2;
		nloop2 = 5;
	}

	hConvergenceCheckOP convergence_checker ( nullptr );
	if ( !option[ basic::options::OptionKeys::abinitio::skip_convergence_check ] ) {
		convergence_checker = hConvergenceCheckOP( new hConvergenceCheck );
	}

	moves::TrialMoverOP trials = trial_large();
	int iteration = 1;
	for ( int lct1 = 1; lct1 <= nloop1; lct1++ ) {
		if ( lct1 > 1 ) trials = trial_small(); //only with short_insert_region!
		for ( int lct2 = 1; lct2 <= nloop2; lct2++, iteration++  ) {
			tr.Debug << "Loop: " << lct1 << "   " << lct2 << std::endl;

			if ( !prepare_loop_in_stage3( pose, iteration, nloop1*nloop2 ) ) return false;

			if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(), "stage_3_iter"+string_of( lct1)+"_"+string_of(lct2),
					false /*fullatom */, true /*fold tree */ ) ) {


				tr.Debug << "  Score stage3 loop iteration " << lct1 << " " << lct2 << std::endl;
				if ( convergence_checker ) {
					moves::TrialMoverOP stage3_trials = stage3_mover( pose, lct1, lct2, trials );
					convergence_checker->set_trials( stage3_trials ); //can be removed late
					moves::WhileMover( stage3_trials, stage3_cycles(), convergence_checker ).apply( pose );
				} else {    //no convergence check -> no WhileMover
					moves::RepeatMover( stage3_mover( pose, lct1, lct2, trials ), stage3_cycles() ).apply( pose );
				}

				if ( numeric::mod( (int)iteration, 2 ) == 0 || iteration > 7 ) recover_low( pose, STAGE_3a );
				recover_low( pose, STAGE_3b );

				get_checkpoints().checkpoint( pose, get_current_tag(), "stage_3_iter"+string_of( lct1)+"_"+string_of(lct2), true /*fold tree */ );
			}//recover_checkpoint
			get_checkpoints().debug( get_current_tag(), "stage_3_iter"+string_of( lct1)+"_"+string_of(lct2), current_scorefxn()( pose ) );

			//   structure_store().push_back( mc_->lowest_score_pose() );
		} // loop 2
	} // loop 1
	return true;
}


// interlaced score2 / score 5 loops
/*! @detail stage4 cycles:
nloop_stage4: iterations
stage4_cycle : trials per  iteration

first iteration: use trial_small()
following iterations: use trial_smooth()
only trial_smooth() if just_smooth_cycles==true

prepare_loop_in_stage4() is called each time before the stage4_cycles_ of trials are started.

start each iteration with the lowest_score_pose. ( mc->recover_low()  in prepare_loop_in_stage4()  )

*/
bool ClassicAbinitio::do_stage4_cycles( pose::Pose &pose ) {
	Size nloop_stage4 = 3;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[corrections::score::cenrot]() ) nloop_stage4=2;

	for ( Size kk = 1; kk <= nloop_stage4; ++kk ) {
		tr.Debug << "prepare ..." << std::endl ;
		if ( !prepare_loop_in_stage4( pose, kk, nloop_stage4 ) ) return false;

		if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(), "stage4_kk_" + ObjexxFCL::string_of(kk), false /* fullatom */, true /* fold_tree */ ) ) {
			moves::TrialMoverOP trials;
			if ( kk == 1 && !just_smooth_cycles_ ) {
				trials = trial_small();
			} else {
				tr.Debug << "switch to smooth moves" << std::endl;
				trials = trial_smooth();
			}

			tr.Debug << "start " << stage4_cycles() << " cycles" << std::endl;
			moves::RepeatMover( stage4_mover( pose, kk, trials ), stage4_cycles() ).apply(pose);
			tr.Debug << "finished" << std::endl;
			recover_low( pose, STAGE_4 );

			get_checkpoints().checkpoint( pose, get_current_tag(), "stage4_kk_" + ObjexxFCL::string_of(kk), true /*fold tree */ );
		}
		get_checkpoints().debug( get_current_tag(), "stage4_kk_" + ObjexxFCL::string_of(kk),  current_scorefxn()( pose ) );

		//don't store last structure since it will be exactly the same as the final structure delivered back via apply
		//  if( kk < nloop_stage4 ) // <-- this line was missing although the comment above was existant.
		//   structure_store().push_back( mc_->lowest_score_pose() );
	}  // loop kk

	if ( option[corrections::score::cenrot] ) {
		//switch to cenrot model
		tr.Debug << "switching to cenrot model ..." << std::endl;
		protocols::simple_moves::SwitchResidueTypeSetMover to_cenrot(chemical::CENTROID_ROT);
		to_cenrot.apply(pose);

		//init pose
		(*score_stage4rot_)( pose );
		pack_rotamers_->score_function(score_stage4rot_sc_);
		pack_rotamers_->apply(pose);

		mc_->reset(pose);
		replace_scorefxn( pose, STAGE_4rot, 0 );
		//mc_->set_temperature(1.0);
		//mc_->set_autotemp(true, 1.0);

		//debug
		//tr.Debug << "starting_energy: " << (*score_stage4rot_)( pose ) << std::endl;
		//tr.Debug << "starting_temperature: " << mc_->temperature() << std::endl;

		for ( Size rloop=1; rloop<=3; rloop++ ) {
			//change vdw weight
			switch (rloop) {
			case 1 :
				score_stage4rot_->set_weight(core::scoring::vdw, score_stage4rot_->get_weight(core::scoring::vdw)/9.0);
				break;
			case 2 :
				score_stage4rot_->set_weight(core::scoring::vdw, score_stage4rot_->get_weight(core::scoring::vdw)*3.0);
				break;
			case 3 :
				score_stage4rot_->set_weight(core::scoring::vdw, score_stage4rot_->get_weight(core::scoring::vdw)*3.0);
				break;
			}

			//stage4rot
			//for (Size iii=1; iii<=100; iii++){
			//pose::Pose startP = pose;
			//tr << "temperature: " << mc_->temperature() << std::endl;
			moves::RepeatMover( stage4rot_mover( pose, rloop, trial_smooth() ), stage4_cycles()/100 ).apply(pose);
			//tr << "delta_rms: " << core::scoring::CA_rmsd( startP, pose ) << std::endl;
			//}
		}
	}

	return true;
}

bool ClassicAbinitio::do_stage5_cycles( pose::Pose &pose ) {//vats

	Size nmoves = 1;
	core::kinematics::MoveMapOP mm_temp( new core::kinematics::MoveMap( *movemap() ) );
	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( mm_temp, temperature_, nmoves) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 2.0 );
	small_mover->angle_max( 'L', 5.0 );

	moves::TrialMoverOP trials( new moves::TrialMover( small_mover, mc_ptr() ) );
	moves::RepeatMover( stage5_mover( pose, trials ), stage5_cycles() ).apply( pose );

	// moves::MoverOP trial( stage5_mover( pose, small_mover ) );
	// Size j;
	// for( j = 1; j <= stage5_cycles(); ++j ) {
	//  trial->apply( pose );
	// }
	mc().reset( pose );
	return true;

}


moves::TrialMoverOP
ClassicAbinitio::stage1_mover( pose::Pose &, moves::TrialMoverOP trials ) {
	return trials;
}


moves::TrialMoverOP
ClassicAbinitio::stage2_mover( pose::Pose &, moves::TrialMoverOP trials ) {
	return trials;
}


moves::TrialMoverOP
ClassicAbinitio::stage3_mover( pose::Pose &, int, int, moves::TrialMoverOP trials ) {
	return trials;
}

moves::TrialMoverOP
ClassicAbinitio::stage4_mover( pose::Pose &, int, moves::TrialMoverOP trials ) {
	return trials;
}

moves::TrialMoverOP
ClassicAbinitio::stage4rot_mover( pose::Pose &, int, moves::TrialMoverOP trials ) {
	if ( trials == trial_small_ ) {
		return trial_small_pack_;
	} else {
		return smooth_trial_small_pack_;
	}
}

moves::TrialMoverOP //vats
ClassicAbinitio::stage5_mover( pose::Pose &, moves::TrialMoverOP trials ) {
	return trials;
}

void ClassicAbinitio::recover_low( core::pose::Pose& pose, StageID stage ){
	if ( contains_stageid( recover_low_stages_, stage ) ) {
		mc_->recover_low( pose );
	}
}

// anything you want to have done before the stages ?
void ClassicAbinitio::replace_scorefxn( core::pose::Pose& pose, StageID stage, core::Real /*intra_stage_progress */ ) {
	// must assume that the current pose is the one to be accepted into the next stage! (this change was necessary for
	// checkpointing to work correctly.

	//intra_stage_progress = intra_stage_progress;
	if ( score_stage1_  && ( stage == STAGE_1 ) ) current_scorefxn( *score_stage1_ );
	if ( score_stage2_  && ( stage == STAGE_2 ) ) current_scorefxn( *score_stage2_ );
	if ( score_stage3a_ && ( stage == STAGE_3a) ) current_scorefxn( *score_stage3a_ );
	if ( score_stage3b_ && ( stage == STAGE_3b) ) current_scorefxn( *score_stage3b_ );
	if ( score_stage4_  && ( stage == STAGE_4 ) ) current_scorefxn( *score_stage4_ );
	if ( score_stage4rot_  && ( stage == STAGE_4rot ) ) current_scorefxn( *score_stage4rot_ );
	if ( score_stage5_  && ( stage == STAGE_5 ) ) current_scorefxn( *score_stage5_ );//vats
	Real temperature( temperature_ );
	if ( stage == STAGE_5 ) temperature = 0.5;
	mc_->set_autotemp( true, temperature );
	mc_->set_temperature( temperature ); // temperature might have changed due to autotemp..
	mc_->reset( pose );
}


moves::TrialMoverOP ClassicAbinitio::trial_large() {
	return ( apply_large_frags_ ? trial_large_ : trial_small_ );
}

moves::TrialMoverOP ClassicAbinitio::trial_small() {
	return trial_small_;
}

moves::TrialMoverOP ClassicAbinitio::trial_smooth() {
	return smooth_trial_small_;
}

// prepare stage1 sampling
bool ClassicAbinitio::prepare_stage1( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_1, 0.5 );
	mc_->set_autotemp( false, temperature_ );
	// mc_->set_temperature( temperature_ ); already done in replace_scorefxn
	// mc_->reset( pose );
	(*score_stage1_)( pose );
	/// Now handled automatically.  score_stage1_->accumulate_residue_total_energies( pose ); // fix this
	return true;
}

bool ClassicAbinitio::prepare_stage2( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_2, 0.5 );

	(*score_stage2_)(pose);
	/// Now handled automatically.  score_stage2_->accumulate_residue_total_energies( pose );
	return true;
}


bool ClassicAbinitio::prepare_stage3( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_3a, 0 );
	//score for this stage is changed in the do_stage3_cycles explicitly
	if ( option[ templates::change_movemap ].user() && option[ templates::change_movemap ] == 3 ) {
		kinematics::MoveMapOP new_mm( new kinematics::MoveMap( *movemap() ) );
		new_mm->set_bb( true );
		set_movemap( new_mm ); // --> store it in movemap_ --> original will be reinstated at end of apply()
	}
	return true;
}


bool ClassicAbinitio::prepare_stage4( core::pose::Pose &pose ) {
	replace_scorefxn( pose, STAGE_4, 0 );
	(*score_stage4_)( pose );
	/// Now handled automatically.  score_stage4_->accumulate_residue_total_energies( pose ); // fix this

	if ( option[ templates::change_movemap ].user() && option[ templates::change_movemap ] == 4 ) {
		kinematics::MoveMapOP new_mm( new kinematics::MoveMap( *movemap() ) );
		new_mm->set_bb( true );
		tr.Debug << "option: templates::change_movemap ACTIVE: set_movemap" << std::endl;
		set_movemap( new_mm ); // --> store it in movemap_ --> original will be reinstated at end of apply()
	}
	return true;
}

bool ClassicAbinitio::prepare_stage5( core::pose::Pose &pose ) {//vats
	// temperature_ = 0.5; //this has to be reset to original temperature!!!
	// no special if-statement in replace_scorefxn...OL
	replace_scorefxn( pose, STAGE_5, 0 );
	(*score_stage5_)( pose );
	return true;
}


bool ClassicAbinitio::prepare_loop_in_stage3( core::pose::Pose &pose/*pose*/, Size iteration, Size total ){
	// interlace score2/score5

	Real chbrk_weight_stage_3a = 0;
	Real chbrk_weight_stage_3b = 0;

	if ( numeric::mod( (int)iteration, 2 ) == 0 || iteration > 7 ) {
		Real progress( iteration );
		chbrk_weight_stage_3a = 0.25 * progress;
		tr.Debug << "select score_stage3a..." << std::endl;
		recover_low( pose, STAGE_3a );
		replace_scorefxn( pose, STAGE_3a, 1.0* iteration/total );
	} else {
		Real progress( iteration );
		chbrk_weight_stage_3b = 0.05 * progress;
		tr.Debug << "select score_stage3b..." << std::endl;
		recover_low( pose, STAGE_3b );
		replace_scorefxn( pose, STAGE_3b, 1.0* iteration/total );
	}

	if ( close_chbrk_ ) {

		set_score_weight( scoring::linear_chainbreak, chbrk_weight_stage_3a , STAGE_3a );
		set_score_weight( scoring::linear_chainbreak, chbrk_weight_stage_3b , STAGE_3b );

	}


	return true;
}

bool ClassicAbinitio::prepare_loop_in_stage4( core::pose::Pose &pose, Size iteration, Size total ){
	replace_scorefxn( pose, STAGE_4, 1.0* iteration/total );

	Real chbrk_weight_stage_4 (iteration*0.5+2.5);

	if ( close_chbrk_ ) {
		set_current_weight( scoring::linear_chainbreak, chbrk_weight_stage_4 );
	}

	return true;
}

//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
moves::MonteCarloOP
ClassicAbinitio::mc_ptr() {
	return mc_;
}


void ClassicAbinitio::output_debug_structure( core::pose::Pose & pose, std::string prefix ) {
	using namespace core::io::silent;

	mc().score_function()( pose );
	Parent::output_debug_structure( pose, prefix );

	if ( option[ basic::options::OptionKeys::abinitio::explicit_pdb_debug ]() ) {
		pose.dump_pdb( prefix + get_current_tag() + ".pdb" );
	}

	if ( option[ basic::options::OptionKeys::abinitio::log_frags ].user() ) {
		std::string filename = prefix + "_" + get_current_tag() + "_" + std::string( option[ basic::options::OptionKeys::abinitio::log_frags ]() );
		utility::io::ozstream output( filename );
		auto& log_frag = dynamic_cast< simple_moves::LoggedFragmentMover& > (*brute_move_large_);
		log_frag.show( output );
		log_frag.clear();
	}

} // ClassicAbinitio::output_debug_structure( core::pose::Pose & pose, std::string prefix )

} //abinitio
} //protocols
