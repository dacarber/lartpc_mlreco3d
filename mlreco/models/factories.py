import torch

def model_dict():

    from . import uresnet
    from . import uresnet_lonely
    from . import uresnet_ppn_chain

    from . import uresnet_clustering
    from . import clustercnn_single
    from . import clustercnn_se

    from . import clusternet
    from . import clustercnn_neural_dbscan
    from . import clustercnn_se
    from . import sparse_occuseg
    from . import sparseoccuseg_gnn

    from . import grappa
    from . import flashmatching_model
    from . import hierarchy
    from . import particle_types

    from . import mink_uresnet
    from . import mink_uresnet_ppn_chain
    from . import mink_singlep
    from . import mink_spice
    from . import vae
    from . import pointnet_gen

    from . import ghost_chain_2
    from . import particle_types
    from . import full_cnn
    from . import full_chain

    # Make some models available (not all of them, e.g. PPN is not standalone)
    models = {
        # Using SCN built-in UResNet
        "uresnet": (uresnet.UResNet, uresnet.SegmentationLoss),
        # Using our custom UResNet
        "uresnet_lonely": (uresnet_lonely.UResNet, uresnet_lonely.SegmentationLoss),
        # URESNET MINKOWSKINET
        "mink_uresnet": (mink_uresnet.UResNet_Chain, mink_uresnet.SegmentationLoss),
        'mink_uresnet_ppn_chain': (mink_uresnet_ppn_chain.UResNetPPN, mink_uresnet_ppn_chain.UResNetPPNLoss), 
        "mink_singlep": (mink_singlep.ParticleImageClassifier, mink_singlep.ParticleTypeLoss),
        "mink_vae": (vae.VAE, vae.ReconstructionLoss), 
        "mink_vae_2": (vae.VAE2, vae.ReconstructionLoss),
        "mink_vae_3": (vae.VAE3, vae.ReconstructionLoss),  
        "mink_spice": (mink_spice.MinkSPICE, mink_spice.SPICELoss), 
        "pointnet_gen": (pointnet_gen.VAE, pointnet_gen.ReconstructionLoss),
        # Chain test for track clustering (w/ DBSCAN)
        #"chain_track_clustering": (chain_track_clustering.Chain, chain_track_clustering.ChainLoss),
        # Chain UResNet and PPN
        "uresnet_ppn_chain": (uresnet_ppn_chain.Chain, uresnet_ppn_chain.ChainLoss),
        # Clustering
        "uresnet_clustering": (uresnet_clustering.UResNet, uresnet_clustering.SegmentationLoss),
        # ClusterUNet Single
        "clustercnn_single": (clustercnn_single.ClusterCNN, clustercnn_single.ClusteringLoss),
        # Colossal ClusterNet Model to Wrap them all
        "clusternet": (clusternet.ClusterCNN, clusternet.ClusteringLoss),
        # Density Loss
        "clustercnn_density": (clustercnn_neural_dbscan.ClusterCNN, clustercnn_neural_dbscan.ClusteringLoss),
        # Spatial Embeddings
        "spatial_embeddings": (clustercnn_se.ClusterCNN, clustercnn_se.ClusteringLoss),
        # OccuSeg
        "occuseg": (sparse_occuseg.SparseOccuSeg, sparse_occuseg.SparseOccuSegLoss),
        # OccuSeg with GNN
        "occuseg_gnn": (sparseoccuseg_gnn.SparseOccuSegGNN, sparseoccuseg_gnn.SparseOccuSegGNNLoss),
        # Spatial Embeddings Lite
        "spatial_embeddings_lite": (clustercnn_se.ClusterCNN2, clustercnn_se.ClusteringLoss),
        # Spatial Embeddings Lovasz free
        "spatial_embeddings_free": (clustercnn_se.ClusterCNN, clustercnn_se.ClusteringLoss),
        # Graph neural network Particle Aggregation (GrapPA)
        "grappa": (grappa.GNN, grappa.GNNLoss),
        # Flashmatching using encoder and gnn
        "flashmatching": (flashmatching_model.FlashMatchingModel, torch.nn.CrossEntropyLoss(reduction='mean')),
        # Particle flow reconstruction with GrapPA (TODO: should be merged with GrapPA)
        "hierarchy_gnn": (hierarchy.ParticleFlowModel, hierarchy.ChainLoss),
        # Particle image classifier
        "particle_type": (particle_types.ParticleImageClassifier, particle_types.ParticleTypeLoss),
        # Deghosting models
        "ghost_chain": (ghost_chain_2.GhostChain2, ghost_chain_2.GhostChain2Loss),
        # CNN chain with UResNet+PPN+SPICE with a single encoder
        "full_cnn": (full_cnn.FullChain, full_cnn.FullChainLoss),
        # Full reconstruction chain, including an option for deghosting
        "full_chain": (full_chain.FullChain, full_chain.FullChainLoss),
    }
    return models


def construct(name):
    models = model_dict()
    if name not in models:
        raise Exception("Unknown model name provided: %s" % name)
    return models[name]
