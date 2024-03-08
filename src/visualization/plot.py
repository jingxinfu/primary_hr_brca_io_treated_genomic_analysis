import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from ..settings  import COLOR_PAlETTE

def responsePlot_highlight_er_low(data,features,ncols):
    df = data[features+['Treatment_Arm','BestResponse','er_status']].copy()
    fig,axs = plt.subplots(2,ncols,figsize=(4*ncols,6),sharex=False,sharey=False)
    order=['Chemo->ICI','ICI->Chemo']
    hue='BestResponse'
    hue_order = ['0-I','II-III']
    for i,y in enumerate(features):
        ax =axs[0,i]
        x = 'BestResponse'
        sns.boxplot(data=df,x=x,y=y,color='white',ax=ax,order=hue_order)
        sns.stripplot(data = df.loc[df.er_status!='Weakly Positive (1-10% cell staining)'],x=x,y=y,hue=hue,ax=ax,order=hue_order,
                      palette=COLOR_PAlETTE[hue],size=7)
    
        sns.stripplot(data = df.loc[df.er_status=='Weakly Positive (1-10% cell staining)'],x=x,y=y,hue=hue,ax=ax,order=hue_order,
                      palette=COLOR_PAlETTE[hue],marker='X', linewidth=1,edgecolor='red',label='ER low',size=7)
        annot = Annotator(ax, pairs=[('0-I','II-III')],data=df, x=x, y=y,order=hue_order)
        annot.configure(test='Mann-Whitney', text_format='full', loc='outside', verbose=0)
        annot._pvalue_format.pvalue_format_string='{:.2g}'
        annot.apply_test()
        annot.annotate()

        ax.set(ylabel='',xlabel='')
        ax.set_title(y.replace('bulkRNA_','').replace('WES_',''),y=1.3)
        if i > 0:
            ax.legend_.remove()
        else:
           handles, labels = ax.get_legend_handles_labels()
           print(handles)
           by_label = dict(zip(labels, handles))
           ax.legend(by_label.values(), by_label.keys(),loc=(-.5,1.15))

        ax=axs[1,i]
        x='Treatment_Arm'
        sns.stripplot(data = df.loc[df.er_status!='Weakly Positive (1-10% cell staining)'],x=x,y=y,hue=hue,ax=ax,order=order,
                      palette=COLOR_PAlETTE[hue],size=7,dodge=True)
    
        sns.stripplot(data = df.loc[df.er_status=='Weakly Positive (1-10% cell staining)'],x=x,y=y,hue=hue,ax=ax,order=order,
                      palette=COLOR_PAlETTE[hue],marker='X', linewidth=1,edgecolor='red',label='ER low',size=7,dodge=True)
        annot = Annotator(ax, pairs=[
            (('Chemo->ICI','0-I'),('Chemo->ICI','II-III')),
            (('ICI->Chemo','0-I'),('ICI->Chemo','II-III')),
            ],
                          data=df, x=x, y=y, hue=hue,order=order,hue_order=hue_order)
        annot.configure(test='Mann-Whitney', text_format='full', loc='outside', verbose=0)
        annot._pvalue_format.pvalue_format_string='{:.2g}'
        annot.apply_test()
        annot.annotate()
        ax.set(ylabel='',xlabel='')
        ax.tick_params(axis='x', labelrotation = 90)
        ax.legend_.remove()

    plt.subplots_adjust(hspace=.5,wspace=.3)
    return fig

def treatmentOverallPlot(data,features,ncols,timepoints):
    df = data[features+['Timepoint','BestResponse','Treatment_Arm','Patient']].copy()
    fig,axs = plt.subplots(2,ncols,figsize=(4*ncols,8),sharex=False,sharey=False)
    assert len(timepoints)==2, "only support 2 timepoints."
    order=['Chemo->ICI','ICI->Chemo']
    hue='Timepoint'
    hue_order = timepoints
    pts = df['Patient'].value_counts()
    selected_pts = pts[pts>1].index
    df = df.loc[df.Patient.isin(selected_pts),:]
    for i,y in enumerate(features):
        # All Treatment
        x='All'
        df[x] = ','.join(df['Treatment_Arm'].unique().tolist())
        x_name = df[x].unique()[0]

        ax=axs[0,i]
        sns.stripplot(data=df,x=x,y=y,hue=hue,palette=COLOR_PAlETTE[hue],ax=ax,hue_order=hue_order,dodge=True)
        ## add pair line
        for _,row in df.pivot_table(index='Patient',columns='Timepoint',values=y).iterrows():
            y0,y1= row[hue_order]
            ax.plot([-.25, .25], [y0,y1], color='black', ls=':', zorder=0)
            
        annot = Annotator(ax, pairs=[
            ((x_name,timepoints[0]),(x_name,timepoints[1])),
            ],data=df, x=x, y=y,hue=hue,hue_order=hue_order)
        annot.configure(test='Wilcoxon', text_format='simple', loc='outside', verbose=1)
        annot.apply_test()
        annot.annotate()
        ax.set(xlabel='',ylabel='')
        ax.set_title(y.replace('bulkRNA_','').replace('WES_',''),y=1.3)
        if i > 0:
            ax.legend_.remove()
        else:
            ax.legend(loc=(-.7,1.15))
        
        # Seperate Treatment
        ax =axs[1,i]
        x='Treatment_Arm'
        sns.stripplot(data=df,x=x,y=y,hue=hue,palette=COLOR_PAlETTE[hue],ax=ax,order=order,hue_order=hue_order,dodge=True)
        ## add pair line
        for _,row in df.loc[df[x]==order[0],:].pivot_table(index='Patient',columns='Timepoint',values=y).iterrows():
            y0,y1= row[hue_order]
            ax.plot([-.25, .25], [y0,y1], color='black', ls=':', zorder=0)
        
        for _,row in df.loc[df[x]==order[1],:].pivot_table(index='Patient',columns='Timepoint',values=y).iterrows():
            y0,y1= row[hue_order]
            ax.plot([.75, 1.25], [y0,y1], color='black', ls=':', zorder=0)
            
        annot = Annotator(ax, pairs=[
            (('Chemo->ICI',timepoints[0]),('Chemo->ICI',timepoints[1])),
            (('ICI->Chemo',timepoints[0]),('ICI->Chemo',timepoints[1])),
            ],
                          data=df, x=x, y=y, hue=hue,order=order,hue_order=hue_order)
        annot.configure(test='Wilcoxon', text_format='simple', loc='outside', verbose=0)
        annot.apply_test()
        annot.annotate()
        ax.set(xlabel='',ylabel='',title='')
        ax.legend_.remove()

    plt.subplots_adjust(hspace=.5,wspace=.3)
    return fig



def treatmentPlot(data,features,ncols,timepoints):
    df = data[features+['Timepoint','BestResponse','Treatment_Arm','Patient']].copy()
    fig,axs = plt.subplots(3,ncols,figsize=(4*ncols,8),sharex=True,sharey=False)
    assert len(timepoints)==2, "only support 2 timepoints."
    order=['Chemo->ICI','ICI->Chemo']
    hue='Timepoint'
    hue_order = timepoints
    for i,y in enumerate(features):
        ax =axs[0,i]
        x='Treatment_Arm'
        sns.stripplot(data=df,x=x,y=y,hue=hue,palette=COLOR_PAlETTE[hue],ax=ax,order=order,hue_order=hue_order,dodge=True)
        ## add pair line
        for _,row in df.loc[df[x]==order[0],:].pivot_table(index='Patient',columns='Timepoint',values=y).iterrows():
            y0,y1= row[hue_order]
            ax.plot([-.25, .25], [y0,y1], color='black', ls=':', zorder=0)
        
        for _,row in df.loc[df[x]==order[1],:].pivot_table(index='Patient',columns='Timepoint',values=y).iterrows():
            y0,y1= row[hue_order]
            ax.plot([.75, 1.25], [y0,y1], color='black', ls=':', zorder=0)
            
        annot = Annotator(ax, pairs=[
            (('Chemo->ICI',timepoints[0]),('Chemo->ICI',timepoints[1])),
            (('ICI->Chemo',timepoints[0]),('ICI->Chemo',timepoints[1])),
            ],
                          data=df, x=x, y=y, hue=hue,order=order,hue_order=hue_order)
        annot.configure(test='Mann-Whitney', text_format='simple', loc='outside', verbose=0)
        annot.apply_test()
        annot.annotate()
        ax.set(xlabel='',ylabel='RCB=0-III')
        ax.set_title(y.replace('bulkRNA_','').replace('WES_',''),y=1.3)
        
        if i > 0:
            ax.legend_.remove()
        else:
            ax.legend(loc=(-.7,1.15))
        # R
        ax=axs[1,i]
        response = '0-I'
        sub_df = df.loc[df.BestResponse==response,:]
        sns.stripplot(data=sub_df,x=x,y=y,hue=hue,palette=COLOR_PAlETTE[hue],ax=ax,order=order,hue_order=hue_order,dodge=True)
        ## add pair line
        for _,row in sub_df.loc[sub_df[x]==order[0],:].pivot_table(index='Patient',columns='Timepoint',values=y).iterrows():
            y0,y1= row[hue_order]
            ax.plot([-.25, .25], [y0,y1], color='black', ls=':', zorder=0)
        
        for _,row in sub_df.loc[sub_df[x]==order[1],:].pivot_table(index='Patient',columns='Timepoint',values=y).iterrows():
            y0,y1= row[hue_order]
            ax.plot([.75, 1.25], [y0,y1], color='black', ls=':', zorder=0)
            
        annot = Annotator(ax, pairs=[
            (('Chemo->ICI',timepoints[0]),('Chemo->ICI',timepoints[1])),
            (('ICI->Chemo',timepoints[0]),('ICI->Chemo',timepoints[1])),
            ],
                          data=df, x=x, y=y, hue=hue,order=order,hue_order=hue_order)
        annot.configure(test='Mann-Whitney', text_format='simple', loc='outside', verbose=0)
        annot.apply_test()
        annot.annotate()
        ax.set(xlabel='',ylabel=f'RCB={response}')
        ax.legend_.remove()
        
        # NR
        ax=axs[2,i]
        response='II-III'
        sub_df = df.loc[df.BestResponse==response,:]
        sns.stripplot(data=sub_df,x=x,y=y,hue=hue,palette=COLOR_PAlETTE[hue],ax=ax,order=order,hue_order=hue_order,dodge=True)
        ## add pair line
        for _,row in sub_df.loc[sub_df[x]==order[0],:].pivot_table(index='Patient',columns='Timepoint',values=y).iterrows():
            y0,y1= row[hue_order]
            ax.plot([-.25, .25], [y0,y1], color='black', ls=':', zorder=0)
        
        for _,row in sub_df.loc[sub_df[x]==order[1],:].pivot_table(index='Patient',columns='Timepoint',values=y).iterrows():
            y0,y1= row[hue_order]
            ax.plot([.75, 1.25], [y0,y1], color='black', ls=':', zorder=0)
        annot = Annotator(ax, pairs=[
            (('Chemo->ICI',timepoints[0]),('Chemo->ICI',timepoints[1])),
            (('ICI->Chemo',timepoints[0]),('ICI->Chemo',timepoints[1])),
            ],
                          data=df, x=x, y=y, hue=hue,order=order,hue_order=hue_order)
        annot.configure(test='Mann-Whitney', text_format='simple', loc='outside', verbose=0)
        annot.apply_test()
        annot.annotate()
        
        ax.set(xlabel='',ylabel=f'RCB={response}')
        ax.tick_params(axis='x', labelrotation = 90)
        ax.legend_.remove()

    plt.subplots_adjust(hspace=.5,wspace=.3)
    return fig

