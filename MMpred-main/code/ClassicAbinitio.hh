// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/ClassicAbinitio.hh
/// @brief header file for ClassicAbinitio protocol
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange
/// @author James Thompson
/// @author Mike Tyka


#ifndef INCLUDED_protocols_abinitio_ClassicAbinitio_hh
#define INCLUDED_protocols_abinitio_ClassicAbinitio_hh

// Unit Headers


#include <protocols/simple_moves/FragmentMover.fwd.hh>

#include <core/types.hh>

#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

#include <protocols/abinitio/Protocol.hh>
#include <protocols/moves/TrialMover.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>
#include <vector>

//000000000000000000000000
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <sstream>

//======================================================
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <fstream>
#include <sstream>
#include <math.h>
// #include <tuple>
#include <string>
//======================================================

//000000000000000000000000

#include <core/fragment/FragSet.fwd.hh>
#include <utility/vector1.hh>

//000000000000000000000000
using namespace std;
extern core::pose::Pose tempPose;
extern core::pose::Pose native_pose;
//000000000000000000000000

namespace protocols {
namespace abinitio {

/// Move these forward declarations to ClassicAbinitio.fwd.hh
class ClassicAbinitio;
typedef utility::pointer::shared_ptr< ClassicAbinitio > ClassicAbinitioOP;

//@ brief The Classic Abinitio protocol from rosetta++
/*!
@ detail
general usage:
ClassicAbinitio  abinitio;
abinitio.init( pose );
...
while(nstruct) {
abinitio.apply( pose );
}

call ClassicAbinitio::register_options() before core::init::init to add relevant options to the applications help

, with the following
stages, all of which uses a different ScoreFunction based on the cen_std.wts in minirosetta_database:

- Stage 1: large (usually 9mer) randomly selected fragment insertions, only VDW term turned on.
Uses score0.wts_patch and runs for either a maximum of 2000 cycles or until all moveable phi/psi values
have been changed.

- Stage 2: large randomly selected fragment insertions, more score terms turned on. Uses score1.wts_patch
and runs for 2000 cycles.

- Stage 3: uses large randomly selected fragment insertions, although the size of the fragment insertions
is tunable via the set_apply_large_frags( bool ) method. Alternates between score2.wts_patch and score5.wts_patch,
running tunable numbers of 2000-cycle iterations between the two scoring functions.

- Stage 4: uses small (usually 3mer) fragment insertions with the fragment selection based on the Gunn cost for
finding local fragment moves. Runs for 4000-cycles and uses score3.wts_patch.

The class implements the basic abinito approach as known from rosetta++. We tried to set this up, such that
behaviour of the protocol can be changed in many different ways ( see, e.g., FoldConstraints ). To be able to change the
behaviour of the protocol easily the class-apply function and methods called therein (e.g., prepare_XXX() / do_XXX_cycles() ) should
not directly change moves or trials. A reference to the currently used score-function should be obtained by
mc().score_function() ...

Behaviour can be changed in the following ways:

use non-classic FragmentMover  --> eg. not uniformly sampled fragments, but using some weighting
--> large and small moves doesn't have to be 3mers and 9mers... use other movers...
---> or other fragets for the "convenience constructor"
use custom trial classes --> overload update_moves()

change sampling behaviour:
overload prepare_XXX() methods: these are called before the cycling for a certain stage begins
overload do_stageX_cycles() : the actual loops over trial-moves ...

change scoring functions:
overload set_default_scores()
weight-changes effective for all stages: set_score_weight()

*/


class ClassicAbinitio : public Protocol {
	typedef Protocol Parent;
public:
	enum StageID {
		ALL_STAGES = 0,
		STAGE_1,
		STAGE_2,
		STAGE_3a,
		STAGE_3b,
		STAGE_4,
		STAGE_4rot,
		STAGE_5
	};
	/// @brief This constructor does not work -- Fix it before using it.
	// constructor: supply mover classes for Fragment Moves
	ClassicAbinitio(
		simple_moves::FragmentMoverOP brute_move_small,
		simple_moves::FragmentMoverOP brute_move_large,
		simple_moves::FragmentMoverOP smooth_move_small,
		int  /*dummy otherwise the two constructors are ambiguous */
	);

	/// @brief constructor: supply fragsets for large and small fragment moves
	ClassicAbinitio(
		core::fragment::FragSetCOP fragset_small,
		core::fragment::FragSetCOP fragset_large,
		core::kinematics::MoveMapCOP movemap
	);

	/// @brief Explicit copy constructor since this class contains OPs of other classes
	ClassicAbinitio( ClassicAbinitio const & src );

	/// @brief Explicit destructor since this class contains OPs of other classes
	~ClassicAbinitio() override;

	//@brief setup moves, mc-object, scores
	//@details can't call this from constructor; virtual functions don't operate
	//until construction has completed.

	void init( core::pose::Pose const& pose ) override;

	//@brief ClassicAbinitio has virtual functions... use this to obtain a new instance

	moves::MoverOP clone() const override;

	//@brief run protocol on pose
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	//@brief return FramgentMover for smooth_small fragment insertions (i.e., stage4 moves)
	simple_moves::FragmentMoverOP smooth_move_small();

	//@brief return FragmentMover for small fragment insertions ( i.e., stage3/4 moves )
	simple_moves::FragmentMoverOP brute_move_small();

	//@brief return FragmentMover for large fragment insertions (i.e., stage1/2 moves )
	simple_moves::FragmentMoverOP brute_move_large();

	//@brief change the movemap ( is propagated to mover-objects )
	//@detail overload if your extension stores additional moves as member variables
	virtual void set_movemap ( core::kinematics::MoveMapCOP mm );

	//@brief returns current movemap
	core::kinematics::MoveMapCOP movemap();

	//@brief set new instances of FragmentMovers
	void set_moves (
		simple_moves::FragmentMoverOP brute_move_small,
		simple_moves::FragmentMoverOP brute_move_large,
		simple_moves::FragmentMoverOP smooth_move_small
	);

	//@brief set new monte-carlo object
	void set_mc( moves::MonteCarloOP );

	//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	moves::MonteCarlo & mc() {
		return *mc_;
	}

	//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	moves::MonteCarlo const & mc() const {
		return *mc_;
	}

	//@brief override cycle setting ( init() -> sets it according to cmd-line options )
	virtual void set_cycles( core::Real increase_cycles = 1.0 );

	//@brief
	Size total_trials() const {
		return total_trials_;
	}

	/// @brief for debugging, one wants to have access to the native pose.
	// void set_native_pose( core::pose::Pose const & pose );

	//@brief set weight - effective for all scoring functions  stage == -1 --> set weight for all stages
	// mod -1 (ignored) ,  in stage3 mod = 1 --> 3a , mod = 2 --> 3b
	virtual void set_score_weight( core::scoring::ScoreType, core::Real setting, StageID stage = ALL_STAGES );

protected:
	//@brief called to notify about changes regarding movers... new movemap / new instances of FragmentMover
	virtual void update_moves();

	//@brief called by init() --- calls all set_default_XXX methods
	virtual void set_defaults( core::pose::Pose const& pose );

	//@brief read out cmd-line options
	virtual void set_default_options();

	//@brief register cmd-line options in option system ( call before core::init::init )
public:
	static void register_options();
	
	//0000000000000000000000000000000000000000000000000000000000000000000000000000000
	simple_moves::ClassicFragmentMoverOP bms, bml, sms;
	
	virtual void RMA_stage1( core::pose::Pose &pose );
	virtual void RMA_stage2( core::pose::Pose &pose );
	virtual void RMA_stage3( core::pose::Pose &pose );
	virtual void RMA_stage4( core::pose::Pose &pose );
	virtual void RMA_stage5( core::pose::Pose &pose );
	
	
	int lenth_of_sequence;
	int N, G_explore, G_maintain, G_enhance_c, Cluster_num, Crowding_radius_exp, Crowding_radius_mer;

	vector< pair<string,int> > parametersVec;
	vector< vector<core::pose::Pose> > Vec_modal_contact;
	
	int K_Kmediods;
	double temp_mean_num_distance;
	void Map_matrix(vector< vector<double> > Map_matrix_pose);
	void Kmediods(Size, vector<vector<double> > Distance_Matrix);
	vector< vector<int> > best_every_cluster;
	
	vector<core::pose::Pose> All_pose;
	vector< core::pose::Pose > Initial_population_(core::pose::Pose &pose);
	
	void CMscore_and_cluster(vector<core::pose::Pose> &storage_pose);
	double DM_score( core::pose::Pose &pose, core::pose::Pose &tempPose );
	void Last_gen_pose(vector< vector<core::pose::Pose> > Vec_modal);
	void all_gen_pose(vector< vector<core::pose::Pose> > Vec_modal);
	//############################我添加的（final enhance）相关的###################
	//simple_moves::FragmentMoverOP FragAssem_;
	int N1;
	ofstream output1;

	void final_enhence();

	void load_initial_population(string const &route);
	vector<core::pose::Pose> Initial_population;
	vector< vector<double> >dscore_map;
	void lock_position(vector<Size> shib);
	double again_calc();
	Size max_vec_mean_index;
	double max_vec_mean;
	void Last_gen_pose();
	
	double Enhance_Distance_score(core::pose::Pose& pose);
	//############################################################################
	
	void Read_parameters();
	void get_parameters(map<string,int> &parametersMap);
	map<string,int> parametersMap;
	
//=======================================================================================	
	double KT_E, KT_C;
	int num_false;
	double Ave_Cscore;
	double Ave_energy;
	double Average_score_of_Populatin(vector<double> &population_score);
	double Average_score_of_Populatin(double &old_score, double &new_score, double &average_score, const Size& population_size);
	
	int frag_;
	simple_moves::FragmentMoverOP FragAssem_;
	
	
	void Multimodal_explore(vector<core::pose::Pose> &storage_pose, int Gen);
	bool Multimodal_maintain(vector<core::pose::Pose> &storage_pose);
	void Multimodal_enhance_c();
	void Multimodal_explore(vector<core::pose::Pose> &storage_pose, vector< double >& population_Cscore, vector< vector<double> > &DMscore_map);
	void Multimodal_maintain(vector<core::pose::Pose> &population, vector< double >& population_energy, vector< double >& population_Cscore);
	vector<core::pose::Pose> Modal_merging(vector<core::pose::Pose> &population, vector<core::pose::Pose> &top_population, vector<core::pose::Pose> &center_population, 
					      vector< vector<core::pose::Pose> > &vec_top_modal, double &cut_dmscore);
	void Multimodal_enhance_c(vector<core::pose::Pose> &population, vector<double> &population_energy, vector<double> &population_Cscore, int enhance_g);
	
	bool boltzmann_accept( double const &targetEnergy, double const &trialEnergy );
	bool boltzmann_accept( double const &targetEnergy, double const &trialEnergy, double const &KT );
	
	vector<double> calculate_population_energy(vector<core::pose::Pose> &population);
	vector<double> calculate_population_Cscore(vector<core::pose::Pose> &population);
	vector< vector<double> >calculate_population_DMscore(vector<core::pose::Pose> &storage_pose);
	
	static double TrialEnergy_;
	static bool isSimilar(double &energy){
		return ( fabs(energy - TrialEnergy_) < 0.0001 );
	}
//=======================================================================================
	struct DistanceType
	{
		int re1;
		int re2;
		float dist;
		double confidence;
		
		DistanceType(){re1=0; re2=0; dist=0; confidence=0;}
	};

	void Read_distance();
	double Distance_score(core::pose::Pose& pose);
	vector<DistanceType> vec_distance;
//0000000000000000000000000000000000000000000000000000000000000000000000000000000
	
protected:
	//@brief construct default monto-carlo object
	virtual void set_default_mc(
		core::pose::Pose const& pose,
		core::scoring::ScoreFunction const& scorefxn
	);

	//@brief assigns default score-functions to the 4 stages--> score0 / score1 / (score2/score5) / score3
	virtual void set_default_scores();

	//@brief currently used score function ( depends on stage )
	core::scoring::ScoreFunction const& current_scorefxn() const;

	//@brief set current scorefunction
	void current_scorefxn( core::scoring::ScoreFunction const& scorefxn );

	//@brief If appropriate for the stage then recover low mc.
	void recover_low( core::pose::Pose& pose, StageID stage );

	//@brief this is called if the scorefxn is replaced with a new one
	virtual void replace_scorefxn( core::pose::Pose& pose, StageID stage, core::Real intra_stage_progress );

	//@brief set individual weight of current scorefunction
	//@details NOTE: does not change the predefined weights, this information is lost each time the score is set
	//  weight will be overwritten with default ->at beginning of stages, each iteration of loop in stage3
	// we could change this behaviour by having a pointer to the current score and making sure that nobody can
	// change the score inside the mc-object (only const accessor )
	void set_current_weight( core::scoring::ScoreType type, core::Real setting );

	//@brief run cycles for different scoring_stages, return number of steps used
	virtual bool do_stage1_cycles( core::pose::Pose &pose );
	virtual bool do_stage2_cycles( core::pose::Pose &pose );
	virtual bool do_stage3_cycles( core::pose::Pose &pose );
	virtual bool do_stage4_cycles( core::pose::Pose &pose );
	virtual bool do_stage5_cycles( core::pose::Pose &pose );//vats

	//@brief returns true if pose is < 3.0 A rms to last pose sent to this function
	// bool convergence_check( core::pose::Pose const & pose );

	//@brief returns the Mover that is applied inside the stage1 loop
	virtual moves::TrialMoverOP
	stage1_mover( core::pose::Pose &pose,  moves::TrialMoverOP trials_in );


	//@brief returns the Mover that is applied inside the stage3 double loop
	virtual moves::TrialMoverOP
	stage2_mover( core::pose::Pose &pose,  moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage3 double loop
	virtual moves::TrialMoverOP
	stage3_mover( core::pose::Pose & pose, int lct1, int lct2,  moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage4 loop
	virtual moves::TrialMoverOP
	stage4_mover( core::pose::Pose &pose, int kk, moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage4 loop
	virtual moves::TrialMoverOP
	stage4rot_mover( core::pose::Pose &pose, int kk, moves::TrialMoverOP trials_in );

	//@brief returns the Mover that is applied inside the stage5 loop
	virtual moves::TrialMoverOP
	stage5_mover( core::pose::Pose &pose, moves::TrialMoverOP trials_in );//vats

	//@brief called by update_moves() creates the instances of TrialMover with the FragmentMoves
	virtual void set_trials();

	//@brief accessor for instances of TrialMover
	moves::TrialMoverOP trial_large();

	//@brief accessor for instances of TrialMover
	moves::TrialMoverOP trial_small();

	//@brief accessor for instances of TrialMover
	moves::TrialMoverOP trial_smooth();

	// anything you want to have done before the stages ?
	//@brief prepare_stageX is called before do_stageX_cycles... overload to change status/scoring/conformation....
	virtual bool prepare_stage1( core::pose::Pose &pose );
	virtual bool prepare_stage2( core::pose::Pose &pose );
	virtual bool prepare_stage3( core::pose::Pose &pose );
	virtual bool prepare_stage4( core::pose::Pose &pose );
	virtual bool prepare_stage5( core::pose::Pose &pose );//vats

	//@brief called in each iteration of inner loop in stage3 before stage3_cycles_ of trials commence
	virtual bool prepare_loop_in_stage3(
		core::pose::Pose&,
		Size, /* loop_iteration*/
		Size  /* total_iterations */
	);

	//@brief called in each iteration of the loop in stage4 before the stage4_cycles_ of trials commence
	virtual bool prepare_loop_in_stage4(
		core::pose::Pose&,
		Size, /* loop_iteration*/
		Size  /* total_iterations */
	);


public:
	//@brief accessor for boolean flag: just_smooth_cycles
	inline
	bool just_smooth_cycles() const {
		return just_smooth_cycles_;
	}

	//@brief Accessor for number of stage1 cycles
	inline
	Size stage1_cycles() const {
		return stage1_cycles_;
	}

	//@brief Accessor for number of stage2
	inline
	Size stage2_cycles() const {
		return stage2_cycles_;
	}

	//@brief Accessor for number of stage3 cycles
	inline
	Size stage3_cycles() const {
		return stage3_cycles_;
	}

	//@brief Setter for number of stage4 cycles
	inline
	void set_stage4_cycles(Size stage4_cycles_new) {
		stage4_cycles_ = stage4_cycles_new;
	}

	//@brief Accessor for number of stage4 cycles
	inline
	Size stage4_cycles() const {
		return stage4_cycles_;
	}

	//@brief Accessor for number of stage5 cycles  //vats
	inline
	Size stage5_cycles() const {
		return stage5_cycles_;
	}


	//@brief query this flag if you have time-intensive stuff and cut it short
	bool bQuickTest() const {
		return bQuickTest_;
	}

	void set_skip_stage1 ( bool setting ) {
		bSkipStage1_ = setting;
	}

	void set_skip_stage2 ( bool setting ) {
		bSkipStage2_ = setting;
	}

protected:
	void output_debug_structure( core::pose::Pose & pose, std::string prefix ) override;

	//@brief obtain currently used monte-carlo object --> use to obtain current score-func: mc().score_function()
	moves::MonteCarloOP mc_ptr();

protected:
	//@brief  cycles_numbers  -- modified by set_cycles(), set_default_options()
	Size stage1_cycles_; //score0
	Size stage2_cycles_; //score1
	Size stage3_cycles_; //score25
	Size stage4_cycles_; //score3
	Size stage5_cycles_; //score3

	core::Real stage4_cycles_pack_rate_;

public:


private:
	//@brief scoring functions -- modified by set_default_scores() / set_score_weight()
	core::scoring::ScoreFunctionOP score_stage1_;  //score0
	core::scoring::ScoreFunctionOP score_stage2_;  //score1
	core::scoring::ScoreFunctionOP score_stage3a_; //score2
	core::scoring::ScoreFunctionOP score_stage3b_; //score5
	core::scoring::ScoreFunctionOP score_stage4_;  //score3
	core::scoring::ScoreFunctionOP score_stage4rot_;  //score_cenrot
	core::scoring::ScoreFunctionOP score_stage4rot_sc_;  //score_cenrot
	core::scoring::ScoreFunctionOP score_stage5_; //score3 //vats lets try score3 first

	//@brief flags
	bool apply_large_frags_;   //above contig_cut2
	bool short_insert_region_; //below contig_cut3
	bool just_smooth_cycles_;
	bool bQuickTest_;  // this flag might land in base-class
	bool close_chbrk_; // closing chainbreaks in old style pose_abinitio

	//@brief a temperature
	core::Real temperature_;

	//@brief movemap --> which dofs can be moved during abinitio
	core::kinematics::MoveMapCOP movemap_;

	//@brief a MonteCarlo object -- set_default_mc() , access: mc()
	moves::MonteCarloOP mc_;

	// Large and small Fragments
	simple_moves::FragmentMoverOP brute_move_small_;
	simple_moves::FragmentMoverOP brute_move_large_;
	simple_moves::FragmentMoverOP smooth_move_small_;

	// TrialMovers
	moves::TrialMoverOP trial_large_;
	moves::TrialMoverOP trial_small_;
	moves::TrialMoverOP smooth_trial_small_;

	moves::TrialMoverOP trial_small_pack_;
	moves::TrialMoverOP smooth_trial_small_pack_;

	// Cenrot Sidechain Mover
	minimization_packing::PackRotamersMoverOP pack_rotamers_;

	Size total_trials_;

public:
	bool bSkipStage1_;
	bool bSkipStage2_;
	bool bSkipStage3_;
	bool bSkipStage4_;
	bool bSkipStage5_;

	utility::vector1< StageID > recover_low_stages_;

private:
	/// @brief Private, unimplemented assignment operator to prevent assignment of this class.
	/// Copy-constructor copying only.
	ClassicAbinitio const & operator = ( ClassicAbinitio const & src );


};

} // abinitio
} // protocols

#endif
