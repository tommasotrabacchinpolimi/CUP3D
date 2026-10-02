//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Fish.h"
#include "StefanFish.h"
#include "FishLibrary.h"
#include "FishShapes.h"
#include "ObstacleVector.h"

#include <Cubism/ArgumentParser.h>

#include <array>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <queue>

CubismUP_3D_NAMESPACE_BEGIN

/**
 * A self-propelled fish whose body is an immersed obstacle on the Cartesian
 * grid, and whose swimming gait is a traveling curvature wave along the spine.
 *
 * Unlike CarlingFish, we do not prescribe the lateral displacement y(s,t)
 * directly. We prescribe curvature kappa(s,t), integrate it with Frenet
 * frames, and rasterize that midline into the characteristic function chi
 * and the deformation velocity udef. Penalization then forces the fluid
 * inside the body to follow U + omega x r + udef. Translation and rotation
 * of the whole fish still come from the hydrodynamics, not from a prescribed
 * path.
 *
 * You can steer the gait without touching the rigid-body solver: Turn(b)
 * sends a bend down the body (yaw from swimming), and a period action
 * changes how fast the tail beats (usually more thrust if T gets shorter).
 * Optional PID flags try to hold the fish near its start pose by scaling
 * amplitude, adding a uniform bend, or wrapping the body to pitch. Do not
 * turn those on together with school control or RL Turn actions — they
 * overwrite the same curvature and fight each other.
 */
class StefanFish: public Fish
{
public:
  /**
   * Reads the usual Fish arguments, then the Stefan-specific ones: tail-beat
   * period -T, phase -phi, amplitudeFactor, height/width profiles, and the
   * CorrectPosition / CorrectPositionZ / CorrectRoll switches. Builds a
   * CurvatureDefinedFishData midline. PID only works if the initial
   * quaternion is identity.
   */
  StefanFish(SimulationData& s, cubism::ArgumentParser& p);

  /**
   * If true, each create() step adjusts alpha (how hard the fish beats,
   * used to hold x) and beta (a uniform extra bend, used to hold y and yaw)
   * so the fish stays near origC in the horizontal plane.
   */
  bool bCorrectPosition;
  /**
   * If true, create() adjusts gamma so the midline wraps onto a cylinder
   * and the fish pitches toward origC in z.
   */
  bool bCorrectPositionZ;
  /**
   * If true, after Obstacle::computeVelocities we strip the roll rate about
   * the long axis and add a small correction that drives roll back to zero.
   * Assumes the default launch pose (quaternion 1,0,0,0).
   */
  bool bCorrectRoll;
  /**
   * Pose the PID treats as "home". Copied from the constructor position;
   * not updated as the fish swims.
   */
  Real origC[3];
  /** How strongly CorrectPosition reacts to y / yaw error. */
  Real wyp;
  /** How strongly CorrectPositionZ reacts to z / pitch error. */
  Real wzp;

  /**
   * Recent estimates of the lab-frame roll axis (xyz + dt). computeVelocities
   * averages a few seconds of this so the roll damper does not chatter.
   */
  std::deque<std::array<Real,4>> r_axis;

  /**
   * Called every step before the Poisson solve, when we need a fresh shape
   * on the grid. If a PID flag is on, we first look at the current pose,
   * write new alpha / beta / gamma, then Fish::create() turns the midline
   * into chi and udef. If the PIDs are off, this is just the usual fish
   * rasterization.
   */
  void create() override;

  /**
   * Lets the base class pick U and omega from fluid momentum inside chi,
   * then — only if bCorrectRoll — removes spin about the body axis and
   * damps remaining roll. That second part is a kinematic override; it is
   * not a hydrodynamic torque.
   */
  virtual void computeVelocities() override;

  /**
   * Dumps the four midline schedulers into Schedulers_*_*.restart and the
   * PID / RL scalars into the shared restart FILE so a later run can
   * continue the same wave, not start from a straight fish.
   */
  virtual void saveRestart(FILE * f) override;

  /** Inverse of saveRestart. Aborts if the scalar block is truncated. */
  virtual void loadRestart(FILE * f) override;

  /**
   * Hands a decision to the midline. lTact is when the wave event should
   * take effect (usually sim.time + half a period). How many numbers you
   * pass chooses the branch inside execute():
   *
   *   one value  — Turn only: a bend b that travels head to tail
   *   three      — Turn plus a period change T_next = T * (1 + a[1]);
   *                the third entry is ignored
   *   five       — the same, plus three torsion samples (needs control_torsion)
   *
   * Two values do nothing, because execute() has no size-2 branch. Turn and
   * period are wave events: call this on the half-period clock, not every
   * CFD step, or the bend never actually travels.
   *
   * There is a sharp edge: if the fish is forced in z (planar runs), we
   * zero a[1] here under the old name "no pitching". That also kills a
   * period action. School control therefore calls action_curvature and
   * action_period directly instead of going through act().
   */
  void act(const Real lTact, const std::vector<Real>& a) const;

  /**
   * Builds the 25-number vector an RL agent sees: lab position and
   * quaternion, wave phase, velocity and spin scaled by T and L, the last
   * two curvature actions, and viscous traction at three skin sensors
   * (nose and both sides of the head).
   */
  std::vector<Real> state() const;

  /**
   * Returns the argument of the traveling sine, wrapped into [0, 2 pi).
   * Useful so an agent can act at a consistent point in the beat.
   */
  Real getPhase(const Real time) const;

  /**
   * Returns the period we are ramping toward (next_period), not necessarily
   * the period right now. RL and school control wait 0.5 * this between
   * actions so Turn has time to propagate.
   */
  Real getLearnTPeriod() const;

  /**
   * Walks the local velocity blocks and returns the index of the one that
   * contains pos. Returns -1 if that point lives on another rank.
   */
  ssize_t holdingBlockID(const std::array<Real,3> pos) const;

  /**
   * Looks up the surface point closest to pSurf on this rank, reads the
   * viscous traction stored there, and MPI-sums so every rank sees the
   * same three components (ranks that do not own the point contribute 0).
   */
  std::array<Real, 3> getShear(const std::array<Real,3> pSurf) const;
};


/**
 * Owns the spine of a StefanFish: curvature (and optional torsion) as
 * functions of arc length, then the 3D midline that Fish rasterizes.
 *
 * The default beat is a Carling-like traveling wave in curvature:
 *
 *   kappa(s,t) = alpha * amplitudeFactor * A(s) * (sin(arg) + rB + beta)
 *
 * A(s) is a smooth envelope, small at the head and large at the tail.
 * sin(arg) is the left-right undulation. rB comes from Turn() and rides
 * down the body as a traveling extra bend. beta is a uniform extra bend
 * used by the planar PID. Frenet integration turns that kappa into
 * (x,y,z) and the deformation velocity. gamma, if set, then wraps the
 * whole planar midline onto a cylinder so the fish pitches.
 */
class CurvatureDefinedFishData : public FishMidlineData
{
 public:
  Real lastTact = 0;
  Real lastCurv = 0;
  Real oldrCurv = 0;

  /** Tail-beat period currently used inside arg(s,t). */
  Real periodPIDval = Tperiod;
  /** How fast that period is changing while we ramp to next_period. */
  Real periodPIDdif = 0;
  bool TperiodPID = false;
  Real lastTime = 0;

  /**
   * time0 and timeshift keep the wave argument continuous when T changes.
   * Without them a sudden period jump would reset the tail phase.
   */
  Real time0 = 0;
  Real timeshift = 0;

  /**
   * Brings the curvature envelope A(s) up from zero over the first period
   * so the fish does not start with a full-amplitude C-shape.
   */
  Schedulers::ParameterSchedulerVector<6>    curvatureScheduler;
  /**
   * Stores Turn decisions as nodes of a wave. gimmeValues() evaluates that
   * wave in a traveling coordinate so a bend injected at the head reaches
   * the tail about one period later.
   */
  Schedulers::ParameterSchedulerLearnWave<7> rlBendingScheduler;

  /**
   * When true, pitching is done with a torsion spline (execute size 5)
   * instead of wrapping the midline with gamma.
   */
  bool control_torsion{false};

  Schedulers::ParameterSchedulerVector<3>    torsionScheduler;
  std::array<Real,3> torsionValues          = {0,0,0};
  std::array<Real,3> torsionValues_previous = {0,0,0};
  Real Ttorsion_start = 0.0;

  /**
   * alpha scales the whole beat (PID uses it to hold streamwise position).
   * beta adds the same extra curvature everywhere (PID yaw / y).
   * gamma is 1/R for the pitch cylinder; leave it 0 for a flat swim.
   */
  Real  alpha     = 1;
  Real dalpha     = 0;
  Real  beta      = 0;
  Real dbeta      = 0;
  Real  gamma     = 0;
  Real dgamma     = 0;

  Schedulers::ParameterSchedulerScalar periodScheduler;
  Real current_period    = Tperiod;
  Real next_period       = Tperiod;
  Real transition_start  = 0.0;
  Real transition_duration = 0.1*Tperiod;

 protected:
  Real * const rK; ///< Curvature along the spine.
  Real * const vK; ///< Time derivative of curvature (needed for udef).
  Real * const rC; ///< Envelope A(s) interpolated from curvatureScheduler.
  Real * const vC;
  Real * const rB; ///< Traveling turn bias from rlBendingScheduler.
  Real * const vB;

  Real * const rT; ///< Torsion along the spine (0 unless control_torsion).
  Real * const vT;
  Real * const rC_T;
  Real * const vC_T;
  Real * const rB_T;
  Real * const vB_T;

 public:
  CurvatureDefinedFishData(Real L, Real T, Real phi, Real _h, const Real _ampFac)
  : FishMidlineData(L, T, phi, _h, _ampFac),
    rK(_alloc(Nm)),vK(_alloc(Nm)), rC(_alloc(Nm)),vC(_alloc(Nm)), rB(_alloc(Nm)),vB(_alloc(Nm)),
    rT(_alloc(Nm)),vT(_alloc(Nm)), rC_T(_alloc(Nm)),vC_T(_alloc(Nm)), rB_T(_alloc(Nm)),vB_T(_alloc(Nm))
    {}

  /**
   * Older helper used when a PID writes a new period immediately. We freeze
   * the current wave argument into timeshift, then install T * periodFac
   * so the sine does not jump. Prefer action_period for scheduled ramps.
   */
  void correctTailPeriod(const Real periodFac, const Real periodVel, const Real t, const Real dt)
  {
    assert(periodFac>0 && periodFac<2); // would be crazy

    const Real lastArg = (lastTime-time0)/periodPIDval + timeshift;
    time0 = lastTime;
    timeshift = lastArg;
    // so that new arg is only constant (prev arg) + dt / periodPIDval
    // with the new l_Tp:
    periodPIDval = Tperiod * periodFac;
    periodPIDdif = Tperiod * periodVel;
    lastTime = t;
    TperiodPID = true;
  }

  /**
   * Interprets an action vector from act(). One number turns. Three numbers
   * turn and change period (the third is unused padding so we hit this
   * branch). Five numbers also set torsion. Anything else, including two
   * numbers, is silently ignored — there is no default case.
   */
  void execute(const Real time, const Real l_tnext, const std::vector<Real>& input) override;

  ~CurvatureDefinedFishData() override
  {
    _dealloc(rK); _dealloc(vK); _dealloc(rC);
    _dealloc(vC); _dealloc(rB); _dealloc(vB);
    _dealloc(rT); _dealloc(vT); _dealloc(rC_T);
    _dealloc(vC_T); _dealloc(rB_T); _dealloc(vB_T);
  }

  /**
   * Fills kappa(s) from the envelope, the traveling sine, Turn, and beta,
   * integrates the Frenet frame to get the midline, then calls
   * performPitchingMotion. This is the function that actually makes the
   * fish undulate each time step.
   */
  void computeMidline(const Real time, const Real dt) override;

  /**
   * If gamma is away from zero, maps the planar (x,y) midline onto a
   * cylinder of radius 1/gamma so the body pitches in z. If gamma is ~0
   * we skip the wrap and leave z = 0.
   */
  void performPitchingMotion(const Real time);

  /**
   * Recomputes normal and binormal from the midline after pitching, so
   * the elliptical cross-section still sits correctly on the spine.
   */
  void recomputeNormalVectors();

  /**
   * Injects a turn. `action` is the bend amplitude (sign picks left vs
   * right). It is attached at the head at time l_tnext and then travels
   * toward the tail; calling this every CFD step resets that wave.
   */
  void action_curvature(const Real time, const Real l_tnext, const Real action)
  {
    rlBendingScheduler.Turn(action, l_tnext);
  }

  /**
   * Asks for a new tail-beat period T * (1 + action), starting at l_tnext.
   * Negative action shortens T, which usually makes the fish swim faster.
   * Warns if the old period PID is still enabled, because that PID will
   * overwrite this ramp.
   */
  void action_period(const Real time, const Real l_tnext, const Real action)
  {
    if (TperiodPID) std::cout << "Warning: PID controller should not be used with RL." << std::endl;
    current_period = periodPIDval;
    next_period = Tperiod * (1 + action);
    transition_start = l_tnext;
  }

  /**
   * Copies three torsion samples onto the spine. They interpolate from the
   * previous samples over half a period. Only meaningful if control_torsion
   * is on.
   */
  void action_torsion(const Real time, const Real l_tnext, const Real * action)
  {
    for (int i = 0 ; i < 3 ; i++)
    {
      torsionValues_previous [i] = torsionValues[i];
      torsionValues[i] = action[i];
    }
    Ttorsion_start = time;
  }

  /**
   * Builds those three torsion samples so the midline twists like a pitch
   * of curvature proportional to `action` (~ 1/R), then forwards them to
   * action_torsion.
   */
  void action_torsion_pitching_radius(const Real time, const Real l_tnext, const Real action)
  {
    const Real sq = 1.0/pow(2.0,0.5);
    const Real ar [3] = { action * ( norX[0     ] * 0.000 + norZ[0     ] * 1.000) ,
                          action * ( norX[Nm/2-1] * sq    + norZ[Nm/2-1] * sq   ) ,
                          action * ( norX[Nm-1  ] * 1.000 + norZ[Nm-1  ] * 0.000) };
    action_torsion(time, l_tnext, ar);
  }
};


CubismUP_3D_NAMESPACE_END
