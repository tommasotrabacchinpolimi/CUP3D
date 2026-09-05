#include "Common.h"
#include "../Simulation.h"
#include "../Obstacles/ObstacleVector.h"
#include "../Obstacles/Sphere.h"
#include <Cubism/ArgumentParser.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <memory>
#include <string>
#include <vector>

using namespace pybind11::literals;
namespace py = pybind11;

CubismUP_3D_NAMESPACE_BEGIN

namespace {

struct ArgvStorage
{
  std::vector<std::string> args;
  std::vector<char *> argv;

  explicit ArgvStorage(const std::vector<std::string> &in)
  {
    args = in;
    if (args.empty())
      args.emplace_back("prg");
    argv.reserve(args.size());
    for (auto &s : args)
      argv.push_back(const_cast<char *>(s.data()));
  }
};

static ArgvStorage kwargsToArgv(const py::kwargs &kwargs)
{
  std::vector<std::string> args{"prg"};
  for (auto item : kwargs) {
    args.emplace_back("-" + py::cast<std::string>(py::str(item.first)));
    args.emplace_back(py::cast<std::string>(py::str(item.second)));
  }
  return ArgvStorage(args);
}

static std::shared_ptr<Sphere> createSphere(SimulationData &sim, py::kwargs kwargs)
{
  if (!kwargs.contains("L") && kwargs.contains("radius"))
    kwargs["L"] = 2.0 * py::cast<double>(kwargs["radius"]);
  auto storage = kwargsToArgv(kwargs);
  cubism::ArgumentParser parser((int)storage.argv.size(), storage.argv.data());
  return std::make_shared<Sphere>(sim, parser);
}

}  // namespace

void bindObstacles(py::module &m)
{
  py::class_<Obstacle, std::shared_ptr<Obstacle>>(m, "Obstacle")
    .def_readwrite("v_imposed", &Obstacle::transVel_imposed);

  py::class_<Sphere, Obstacle, std::shared_ptr<Sphere>>(m, "SphereObstacle")
    .def(py::init([](SimulationData &sim, py::kwargs kwargs) {
           return createSphere(sim, std::move(kwargs));
         }),
         "sim"_a)
    .def_readonly("radius", &Sphere::radius)
    .def_readwrite("umax", &Sphere::umax)
    .def_readwrite("tmax", &Sphere::tmax)
    .def_readwrite("accel_decel", &Sphere::accel_decel)
    .def_readwrite("bHemi", &Sphere::bHemi);
}

void pySimulationAddObstacle(Simulation &s, std::shared_ptr<Obstacle> obstacle)
{
  s.sim.obstacle_vector->addObstacle(std::move(obstacle));
}

void pySimulationParseAndAddObstacle(Simulation &S, pybind11::object obstacle_args)
{
  if (py::isinstance<Sphere>(obstacle_args)) {
    S.sim.obstacle_vector->addObstacle(py::cast<std::shared_ptr<Sphere>>(obstacle_args));
    return;
  }
  if (py::isinstance<Obstacle>(obstacle_args)) {
    S.sim.obstacle_vector->addObstacle(py::cast<std::shared_ptr<Obstacle>>(obstacle_args));
    return;
  }
  throw std::invalid_argument(py::str(obstacle_args));
}

CubismUP_3D_NAMESPACE_END
