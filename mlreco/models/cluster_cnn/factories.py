from mlreco.models.cluster_cnn import losses


def backbone_dict():
    """
    returns dictionary of clustering models
    """
    from mlreco.models.layers import uresnet
    from mlreco.models.layers import fpn

    models = {
        # Encoder-Decoder Type Backbone Architecture.
        "uresnet": uresnet.UResNet,
        "fpn": fpn.FPN
    }

    return models


def gs_kernel_dict():
    '''
    Returns dictionary of kernel function models.
    '''
    from . import gs_kernels
    kernels = {
        'default': gs_kernels.DefaultKernel,
        'bilinear': gs_kernels.BilinearKernel,
        'bilinear_mlp': gs_kernels.BilinearNNKernel,
        'mixed': gs_kernels.MixedKernel
    }
    return kernels


def gs_kernel_construct(cfg):
    models = gs_kernel_dict()
    name = cfg['name']
    num_features = cfg['num_features']
    args = cfg.get('args', {})
    if not name in models:
        raise Exception("Unknown kernel function name provided")
    return models[name](num_features=num_features, **args)


def cluster_model_dict():
    '''
    Returns dictionary of implemented clustering layers.
    '''
    from . import spatial_embeddings
    from . import graph_spice
    from mlreco.mink.cluster.embeddings import SPICE as MinkSPICE
    models = {
        "spice_cnn": spatial_embeddings.SpatialEmbeddings1,
        "spice_cnn_me": MinkSPICE,
        "spice_cnn_lite": spatial_embeddings.SpatialEmbeddingsLite,
        "graph_spice_embedder": graph_spice.GraphSPICEEmbedder,
        "graph_spice_geo_embedder": graph_spice.GraphSPICEGeoEmbedder
        # "graphgnn_spice": graphgnn_spice.SparseOccuSegGNN
    }
    return models


def spice_loss_dict():
    '''
    Returns dictionary of various clustering losses with enhancements.
    '''
    from . import losses
    # from .graphgnn_spice import SparseOccuSegGNNLoss
    loss = {
        # Hyperspace Clustering Losses
        'single': losses.single_layers.DiscriminativeLoss,
        #'multi': losses.multi_layers.MultiScaleLoss,
        #'multi-weighted': losses.multi_layers.DistanceEstimationLoss3,
        #'multi-repel': losses.multi_layers.DistanceEstimationLoss2,
        #'multi-distance': losses.multi_layers.DistanceEstimationLoss,
        # SPICE Losses
        'se_bce': losses.spatial_embeddings.MaskBCELoss2,
        'se_bce_ellipse': losses.spatial_embeddings.MaskBCELossBivariate,
        'se_lovasz': losses.spatial_embeddings.MaskLovaszHingeLoss,
        'se_lovasz_inter': losses.spatial_embeddings.MaskLovaszInterLoss,
        'se_focal': losses.spatial_embeddings.MaskFocalLoss,
        'se_multivariate': losses.spatial_embeddings.MultiVariateLovasz,
        'se_ce_lovasz': losses.spatial_embeddings.CELovaszLoss,
        'se_lovasz_inter_2': losses.spatial_embeddings.MaskLovaszInterLoss2,
        'se_lovasz_inter_bc': losses.spatial_embeddings.MaskLovaszInterBC,
        # SPICE Losses Vectorized
        'se_vectorized': losses.spatial_embeddings_fast.SPICELoss,
        'se_vectorized_inter': losses.spatial_embeddings_fast.SPICEInterLoss,
        'graph_spice_edge_loss': losses.gs_embeddings.NodeEdgeHybridLoss,
        'graph_spice_loss': losses.gs_embeddings.GraphSPICEEmbeddingLoss
        # 'graphgnn_spice_loss': SparseOccuSegGNNLoss
    }
    return loss


def backbone_construct(name):
    models = backbone_dict()
    if not name in models:
        raise Exception("Unknown backbone architecture name provided")
    return models[name]


def cluster_model_construct(cfg, name):
    models = cluster_model_dict()
    if not name in models:
        raise Exception("Unknown clustering model name provided")
    return models[name](cfg)


def spice_loss_construct(name):
    loss_fns = spice_loss_dict()
    if not name in loss_fns:
        raise Exception("Unknown clustering loss function name provided")
    return loss_fns[name]
