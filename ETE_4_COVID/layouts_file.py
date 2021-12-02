from ete4 import NodeStyle
from ete4.smartview.ete.treestyle import TreeStyle
from ete4.smartview.ete.faces import RectFace, CircleFace, SeqMotifFace, TextFace, PieChartFace
#from ete4 import PieChartFace #NO IMPLEMENTADO PARESE
from collections import Counter
import ast

def get_layout_clade():
    
    col_and_clades={'S':'#eba407','L':'#13a54c','V':'#ddeccb','G':'#142733','GK':'#4c32cf','GH':'#f73630','GR':'#42466a','GV':'#8b2f7d','GRY':'#3c82db','O':'#906175'}    
    lineages=[]

    def layout_fn(node):
        leaf_style=NodeStyle()

        if node.is_leaf():
            try:
                leaf_style['fgcolor']=col_and_clades[node.props['clade']]
                leaf_style["size"] = 5
                #leaf_style["bgcolor"] =col_and_clades[node.props['clade']]
                node.set_style(leaf_style)
                node.add_face(TextFace(node.props['clade'],  padding_x=15), column=0, position="aligned")
            except Exception:
                #node.set_style(leaf2_style)
                #node.add_face(TextFace("Removed sample from GISAID",  padding_x=15), column=0, position="branch_right")
                pass

        else:

            if 'clade_frequency' in node.props:
                   
                clades=ast.literal_eval(node.props['clade_frequency'])
                clades_dict={c.split('_')[0]: [int(c.split('_')[1]),col_and_clades[c.split('_')[0]]] for c in clades}

                data=[ [key, value[0], value[1]] for key,value in clades_dict.items()]
                pie = PieChartFace(20, data=data)
                position = "branch_right"
                node.add_face(pie, position=position, column=1, collapsed_only=(not node.is_leaf()))

                values=[v[1] for v in data]
                try:
                    major_clade=[i[0] for i in data if max(values)==i[1]]
                    node.add_face(TextFace(major_clade,  padding_x=5),column=2, position="branch_right",collapsed_only=(not node.is_leaf()))

                except Exception:
                    major_clade='Nonexistence'
                    node.add_face(TextFace(major_clade,  padding_x=5),column=2, position="branch_right",collapsed_only=(not node.is_leaf()))

            else:
                leaf_style['fgcolor']='black'
                leaf_style["size"] = 5
                leaf_style['shape'] = 'square'
                node.set_style(leaf_style)
                node.add_face(TextFace("UNKNOWN",  padding_x=15), column=2, position="branch_right",collapsed_only=True)


    layout_fn.__name__ = 'GISAID Clades'
    layout_fn.contains_aligned_face = False
    return layout_fn

def get_layout_lineage():
    
    def layout_fn(node):
        
        if node.is_leaf():
            try:
                node.add_face(TextFace(node.props['pango_lineage'],  padding_x=15), column=1, position="aligned") 
            except Exception:
                node.add_face(TextFace("Removed sample from GISAID",  padding_x=15), column=1, position="aligned")
            
        
        
    layout_fn.__name__ = 'Pango lineages'
    layout_fn.contains_aligned_face = True 
    return layout_fn

def get_layout_continent():
    
    col_and_continents={'Asia':'#ac6af4','Africa':'#be3880','Europe':'#76ab86','North America':'#d49e3c','Oceania':'#4c32cf',
                        'South America':'#f73630','Caribbean':'#42466a','Central America':'#8b2f7d'}

    padding_x = padding_y = 1

    def layout_fn(node):
        
        if node.is_leaf():
            try:
                continent=node.props['location'].split("/")[0].strip()
                #node.add_face(TextFace(continent,  padding_x=15), column=1, position="branch_right")
                node_rect=RectFace(50, 75, col_and_continents[continent], padding_x=padding_x, padding_y=padding_y, text=continent) #, rotate_text = False
                node.add_face(node_rect, column=4, position="aligned")
            except Exception:
                pass

        else:
            #--Internal nodes already have the property.
            if 'region_frequency' in node.props:
                   
                regions=ast.literal_eval(node.props['region_frequency'])
                regions_dict={c.split('_')[0]: [int(c.split('_')[1]),col_and_continents[c.split('_')[0]]] for c in regions}

                data=data=[ [key, value[0], value[1]] for key,value in regions_dict.items()]
                pie = PieChartFace(20, data=data)
                position = "branch_right"
                node.add_face(pie, position=position, column=4, collapsed_only=(not node.is_leaf()))

                values=[v[1] for v in data]
                try:
                    major_clade=[i[0] for i in data if max(values)==i[1]]
                    node.add_face(TextFace(major_clade,  padding_x=5),column=5, position="branch_right",collapsed_only=(not node.is_leaf()))

                except Exception:
                    major_clade='Nonexistence'
                    node.add_face(TextFace(major_clade,  padding_x=5),column=5, position="branch_right",collapsed_only=(not node.is_leaf()))

            else:
                node.add_face(TextFace("UNKNOWN",  padding_x=15), column=5, position="branch_right",collapsed_only=True)
        
        
    layout_fn.__name__ = 'Continents'
    layout_fn.contains_aligned_face = True # declare aligned Faces
    return layout_fn

def get_layout_who_label():
    
    col_and_labels={'Alpha':'#e51ffe','Beta':'#7a77a7','Gamma':'#994a53','Delta':'#984bd9','Kappa':'#5305c4',
                    'Epsilon':'#0cb0d3','Eta':'#0562db','Iota':'#07538e','Lambda':'#510e43','Mu':'#da7499', 'uncorrelated' :'#bbb9bb'}
    

    lineages=[]
    def layout_fn(node):
        leaf_style=NodeStyle()
        if node.is_leaf():
            try:
                leaf_style['fgcolor']=col_and_labels[node.props['who_label']]
                leaf_style["size"] = 5
                #leaf_style["bgcolor"] =col_and_labels[node.props['who_label']]
                node.set_style(leaf_style)
                node.add_face(TextFace(node.props['who_label'],  padding_x=15), column=6, position="aligned")
            except Exception:
                #node.set_style(leaf2_style)
                #node.add_face(TextFace("Removed sample from GISAID",  padding_x=15), column=0, position="branch_right")
                pass

        else:

            if 'who_label_frequency' in node.props:
                  
                who_labels=ast.literal_eval(node.props['who_label_frequency'])
                who_labels_dict={c.split('_')[0]: [int(c.split('_')[1]),col_and_labels[c.split('_')[0]]] for c in who_labels}

                data=[ [key, value[0], value[1]] for key,value in who_labels_dict.items()]
                pie = PieChartFace(20, data=data)
                position = "branch_right"
                node.add_face(pie, position=position, column=7, collapsed_only=(not node.is_leaf()))
                values=[v[1] for v in data]
                try:
                    major_clade=[i[0] for i in data if max(values)==i[1]]
                    node.add_face(TextFace(major_clade,  padding_x=5),column=8, position="branch_right",collapsed_only=(not node.is_leaf()))

                except Exception:
                    major_clade='Nonexistence'
                    node.add_face(TextFace(major_clade,  padding_x=5),column=8, position="branch_right",collapsed_only=(not node.is_leaf()))


            else:
                leaf_style['fgcolor']='black'
                leaf_style["size"] = 5
                leaf_style['shape'] = 'square'
                node.set_style(leaf_style)
                node.add_face(TextFace("UNKNOWN",  padding_x=15), column=8, position="branch_right",collapsed_only=True)


    layout_fn.__name__ = 'WHO Labels'
    layout_fn.contains_aligned_face = True
    return layout_fn


def get_layout_nextstrain_clade():
    
    col_and_labels={'20I(Alpha, V1)':'#55cd81','20H(Beta, V2)':'#a91e24','20J(Gamma, V3)':'#13bf8e','21A/21I/21J':'#a98524','21B(Kappa)':'#7b8da0',
                    '21C(Epsilon)':'#433e7a','21C(Epsilon)':'#10119c','21D(Eta)':'#84b3c9','21F(Iota)':'#c548b9','21G(Lambda)':'#d1c897', '21H(Mu)' :'#0fa973', 
                    '20E' :'#9138c3','20B/S732A' :'#3dfcd3','20A/S126A' :'#0aaab0','20A.EU2' :'#c07d11','20A/S439K' :'#cba475','20A/S98F' :'#bbb9bb',
                    '20C/S80Y' :'#a08451','20B/S626S' :'#9e5484','20B/S1122L' :'#0a6863','19B' :'#cd63e1','19A' :'#8f73b7','uncorrelated':'#bbb9bb'}
    lineages=[]
    def layout_fn(node):
        leaf_style=NodeStyle()
        if node.is_leaf():
            try:
                leaf_style['fgcolor']=col_and_labels[node.props['nextstrain_clade']]
                leaf_style["size"] = 5
                #leaf_style["bgcolor"] =col_and_labels[node.props['who_label']]
                node.set_style(leaf_style)
                node.add_face(TextFace(node.props['nextstrain_clade'],  padding_x=15), column=8, position="aligned")
            except Exception:
                #node.set_style(leaf2_style)
                #node.add_face(TextFace("Removed sample from GISAID",  padding_x=15), column=0, position="branch_right")
                pass

        else:

            if 'nextstrain_frequency' in node.props:
                  
                next_labels=ast.literal_eval(node.props['nextstrain_frequency'])
                next_labels_dict={c.split('_')[0]: [int(c.split('_')[1]),col_and_labels[c.split('_')[0]]] for c in next_labels}
                
                data=[ [key, value[0], value[1]] for key,value in next_labels_dict.items()]
                pie = PieChartFace(20, data=data)
                position = "branch_right"
                node.add_face(pie, position=position, column=10, collapsed_only=(not node.is_leaf()))

                values=[v[1] for v in data]
                try:
                    major_clade=[i[0] for i in data if max(values)==i[1]]
                    node.add_face(TextFace(major_clade,  padding_x=5),column=11, position="branch_right",collapsed_only=(not node.is_leaf()))

                except Exception:
                    major_clade='Nonexistence'
                    node.add_face(TextFace(major_clade,  padding_x=5),column=11, position="branch_right",collapsed_only=(not node.is_leaf()))

                


            else:
                leaf_style['fgcolor']='black'
                leaf_style["size"] = 5
                leaf_style['shape'] = 'square'
                node.set_style(leaf_style)
                node.add_face(TextFace("UNKNOWN",  padding_x=5), column=11, position="branch_right",collapsed_only=True)


    layout_fn.__name__ = 'Nextstrain Clades'
    layout_fn.contains_aligned_face = True
    return layout_fn