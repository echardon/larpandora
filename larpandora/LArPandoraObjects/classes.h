#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"

#include "lardataobj/RecoBase/PFParticle.h"

#include "larpandora/LArPandoraObjects/PFParticleMetadata.h"

namespace {
    struct dictionary {
        art::Assns<larpandoraobj::PFParticleMetadata,recob::PFParticle,void> assnMetaToParticle;
        art::Assns<recob::PFParticle,larpandoraobj::PFParticleMetadata,void> assParticleToMeta;
        art::Assns<recob::PFParticle,bool,void> assParticleToBool;
        art::Assns<bool,recob::PFParticle,void> assBoolToParticle;

        art::Wrapper<larpandoraobj::PFParticleMetadata> wrapperMeta;
        art::Wrapper<std::vector<larpandoraobj::PFParticleMetadata> > wrapperMetaVector;
        art::Wrapper<art::Assns<larpandoraobj::PFParticleMetadata,recob::PFParticle,void> > wrapperAssnMetaToParticle;
        art::Wrapper<art::Assns<recob::PFParticle,larpandoraobj::PFParticleMetadata,void> > wrapperAssnParticleToMeta;
        art::Wrapper<art::Assns<recob::PFParticle,bool,void> > wrapperAssnParticleToBool;
        art::Wrapper<art::Assns<bool,recob::PFParticle,void> > wrapperAssnBoolToParticle;
    };
}
