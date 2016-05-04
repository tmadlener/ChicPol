#ifndef BKGHISTOPRODUCER_FACTORY_H_
#define BKGHISTOPRODUCER_FACTORY_H__

#include "BkgHistoProducer.h"

/** Factory that returns a BkgHistoProducer object. */
class BkgHistoProducerFactory {
public:
  /**< return a (pure virtual) instance of a BkgHistoProducer. */
  IBkgHistoProducer* createBkgHistoProducer(const int state);
};

IBkgHistoProducer* BkgHistoProducerFactory::createBkgHistoProducer(const int state)
{
  StateT stateT = stateTfromInt(state);
  switch (stateT) {
  case Jpsi: return new BkgHistoProducer<Jpsi>();
  case Psi2S: return new BkgHistoProducer<Psi2S>();
  case Chic: return new BkgHistoProducer<Chic>();
  case undefined:
    std::cerr << "undefined state in BkgHistoProducerFactory! Returning nullptr" << std::endl;
  default: return NULL; // nullptr;
  }
}

#endif
