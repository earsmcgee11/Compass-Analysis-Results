
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

# Configure page
st.set_page_config(
    page_title="CD4 CD5 Hi/Lo Metabolic Pathway Explorer",
    page_icon="🧬",
    layout="wide"
)

# Title and description
st.title("🧬 CD4 CD5 Hi/Lo Metabolic Pathway Explorer")
st.markdown("""
**Explore differential metabolic flux between CD5 hi and CD5 lo CD4+ T cells**
- Data from Compass metabolic flux analysis
- **ALL pathways** ordered by effect size and significance
- Gene associations from Mouse-GEM metabolic model
- Blue means higher Cohen's D value in CD5lo, Red means higher Cohen's D value in CD5hi
- Cohen's D measures 'how different' two groups are, so a value of 0 means no difference between CD5hi and CD5lo. 
- A Cohen's D value of +2 would mean CD5hi cells have higher metabolic flux in this reaction/pathway, and -2 would mean CD5lo cells have higher metabolic flux in this reaction/pathway.
""")

@st.cache_data
def load_data():
    """Load the comprehensive pathway data"""
    try:
        # Try the new comprehensive file first
        df = pd.read_csv('all_pathways_comprehensive.csv')
        return df
    except FileNotFoundError:
        try:
            # Fallback to old file
            df = pd.read_csv('top_pathways_comprehensive.csv')
            return df
        except FileNotFoundError:
            st.error("Data file not found. Please upload all_pathways_comprehensive.csv")
            return None

# Load data
data = load_data()

if data is not None:
    # Sidebar filters
    st.sidebar.header("🔍 Filters")
    
    # Direction filter
    direction_options = ['All'] + list(data['pathway_direction'].unique())
    selected_direction = st.sidebar.selectbox("Pathway Direction", direction_options)
    
    # Significance filter
    show_significant_only = st.sidebar.checkbox("Show only significant reactions (p < 0.05)", value=False)
    
    # Gene filter
    show_with_genes_only = st.sidebar.checkbox("Show only reactions with known genes", value=False)
    
    # Effect size threshold
    min_effect_size = st.sidebar.slider("Minimum |Effect Size|", 0.0, 5.0, 0.0, 0.1)
    
    # Apply filters
    filtered_data = data.copy()
    
    if selected_direction != 'All':
        filtered_data = filtered_data[filtered_data['pathway_direction'] == selected_direction]
    
    if show_significant_only:
        filtered_data = filtered_data[filtered_data['significant'] == True]
    
    if show_with_genes_only:
        filtered_data = filtered_data[filtered_data['n_genes'] > 0]
    
    # Main content area
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.header("📋 Pathways")
        
        # Get unique pathways from filtered data with stats
        pathways = filtered_data.groupby('pathway').agg({
            'pathway_direction': 'first',
            'pathway_median_d': 'first',
            'significant': 'sum',
            'reaction_id': 'count'
        }).reset_index()
        
        pathways.columns = ['pathway', 'direction', 'median_d', 'n_significant', 'n_total']
        pathways['pct_significant'] = (pathways['n_significant'] / pathways['n_total'] * 100).round(1)
        
        # Filter by effect size
        pathways = pathways[abs(pathways['median_d']) >= min_effect_size]
        
        # Sort pathways by effect size (absolute value) then by % significant
        pathways['abs_median_d'] = abs(pathways['median_d'])
        pathways = pathways.sort_values(['abs_median_d', 'pct_significant'], 
                                       ascending=[False, False])
        pathways = pathways.drop('abs_median_d', axis=1)
        
        st.write(f"**{len(pathways)} pathways match filters**")
        
        # Pathway selector with enhanced display
        selected_pathways = []
        for _, pathway_row in pathways.iterrows():
            pathway_name = pathway_row['pathway']
            direction_emoji = "🔴" if pathway_row['direction'] == 'CD5 hi' else "🔵"
            
            # Enhanced label with effect size and significance
            effect_size = pathway_row['median_d']
            pct_sig = pathway_row['pct_significant']
            n_sig = pathway_row['n_significant']
            n_total = pathway_row['n_total']
            
            # Truncate pathway name and add stats
            short_name = pathway_name[:35] + "..." if len(pathway_name) > 35 else pathway_name
            checkbox_label = f"{direction_emoji} {short_name} (d={effect_size:+.2f}, {pct_sig:.0f}% sig)"
            
            if st.checkbox(checkbox_label, key=f"pathway_{pathway_name}"):
                selected_pathways.append(pathway_name)
        
        # Show selected pathway summary
        if selected_pathways:
            st.subheader("📊 Selected Pathways")
            for pathway in selected_pathways:
                pathway_info = pathways[pathways['pathway'] == pathway].iloc[0]
                direction_color = "🔴" if pathway_info['direction'] == 'CD5 hi' else "🔵"
                st.write(f"{direction_color} **{pathway[:40]}...**")
                st.write(f"   Effect: {pathway_info['median_d']:+.2f} | Significant: {pathway_info['n_significant']}/{pathway_info['n_total']} ({pathway_info['pct_significant']:.1f}%)")
    
    with col2:
        st.header("🧪 Reactions")
        
        if selected_pathways:
            # Filter data for selected pathways
            pathway_data = filtered_data[filtered_data['pathway'].isin(selected_pathways)]
            
            # Sort by effect size (absolute value)
            pathway_data = pathway_data.sort_values('cohens_d', key=abs, ascending=False)
            
            # Display reactions table
            st.write(f"**{len(pathway_data)} reactions in selected pathway(s)**")
            
            # Create display dataframe
            display_df = pathway_data[['reaction_id', 'pathway', 'cohens_d', 'p_value', 'significant', 'genes', 'n_genes', 'reaction_name']].copy()
            display_df['pathway_short'] = display_df['pathway'].str[:25] + "..."
            display_df['genes_short'] = display_df['genes'].str[:40] + "..."
            display_df['cohens_d'] = display_df['cohens_d'].round(2)
            display_df['p_value'] = display_df['p_value'].apply(lambda x: f"{x:.2e}")
            display_df['significant'] = display_df['significant'].map({True: "✓", False: "✗"})
            
            # Select columns to show
            show_columns = ['reaction_id', 'pathway_short', 'cohens_d', 'p_value', 'significant', 'n_genes', 'genes_short']
            
            # Display table
            st.dataframe(
                display_df[show_columns],
                column_config={
                    "reaction_id": "Reaction ID",
                    "pathway_short": "Pathway",
                    "cohens_d": "Cohen's d",
                    "p_value": "p-value",
                    "significant": "Sig",
                    "n_genes": "# Genes",
                    "genes_short": "Associated Genes"
                },
                hide_index=True,
                use_container_width=True
            )
            
            # Reaction details
            st.subheader("🔍 Reaction Details")
            selected_reaction = st.selectbox(
                "Select a reaction for details:",
                options=pathway_data['reaction_id'].tolist(),
                format_func=lambda x: f"{x} (d={pathway_data[pathway_data['reaction_id']==x]['cohens_d'].iloc[0]:+.2f})"
            )
            
            if selected_reaction:
                reaction_info = pathway_data[pathway_data['reaction_id'] == selected_reaction].iloc[0]
                
                col_a, col_b = st.columns(2)
                
                with col_a:
                    st.write(f"**Reaction:** {selected_reaction}")
                    st.write(f"**Name:** {reaction_info['reaction_name']}")
                    st.write(f"**Pathway:** {reaction_info['pathway']}")
                    st.write(f"**Effect Size:** {reaction_info['cohens_d']:+.2f}")
                    st.write(f"**p-value:** {reaction_info['p_value']:.2e}")
                    st.write(f"**Significant:** {'Yes' if reaction_info['significant'] else 'No'}")
                
                with col_b:
                    st.write(f"**Associated Genes ({reaction_info['n_genes']}):**")
                    if reaction_info['genes'] and reaction_info['genes'] != '':
                        genes_list = reaction_info['genes'].split('; ')
                        for gene in genes_list:
                            st.write(f"- {gene}")
                    else:
                        st.write("No genes associated")
                    
                    if reaction_info['ec_number'] and reaction_info['ec_number'] != 'No EC':
                        st.write(f"**EC Number:** {reaction_info['ec_number']}")
        
        else:
            st.info("👈 Select one or more pathways from the left panel to view reactions")
    
    # Summary statistics
    st.header("📊 Summary Statistics")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Total Pathways", len(data['pathway'].unique()))
    
    with col2:
        st.metric("Total Reactions", len(data))
    
    with col3:
        st.metric("Significant Reactions", len(data[data['significant'] == True]))
    
    with col4:
        st.metric("Reactions with Genes", len(data[data['n_genes'] > 0]))
    
    # Pathway overview plot
    st.header("📈 Pathway Overview")
    
    # Create pathway summary for plotting
    pathway_summary = data.groupby(['pathway', 'pathway_direction']).agg({
        'pathway_median_d': 'first',
        'significant': 'sum',
        'reaction_id': 'count'
    }).reset_index()
    
    pathway_summary.columns = ['pathway', 'direction', 'median_d', 'n_significant', 'n_total']
    pathway_summary['pct_significant'] = pathway_summary['n_significant'] / pathway_summary['n_total'] * 100
    
    # Create interactive plot
    fig = px.scatter(
        pathway_summary,
        x='median_d',
        y='pct_significant',
        color='direction',
        size='n_total',
        hover_data=['pathway', 'n_significant', 'n_total'],
        title="Pathway Effect Size vs Significance (All Pathways)",
        labels={
            'median_d': "Median Cohen's d",
            'pct_significant': "% Significant Reactions",
            'direction': "Direction"
        },
        color_discrete_map={'CD5 hi': 'red', 'CD5 lo': 'blue'}
    )
    
    fig.update_layout(height=500)
    st.plotly_chart(fig, use_container_width=True)

else:
    st.error("Please upload the data file to continue.")
