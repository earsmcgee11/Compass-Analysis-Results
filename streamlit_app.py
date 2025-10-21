import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

# Configure page
st.set_page_config(
    page_title="Metabolic Pathway Explorer",
    page_icon="🧬",
    layout="wide"
)

# Title and description
st.title("Metabolic Pathway Explorer")

# Add cache clear button in sidebar
if st.sidebar.button("🔄 Clear Cache"):
    st.cache_data.clear()
    st.sidebar.success("Cache cleared! Please refresh the page.")


@st.cache_data
def load_cd4_data():
    """Load CD4 comprehensive pathway data"""
    try:
        df = pd.read_csv('all_pathways_comprehensive.csv')
        st.sidebar.success(f"✅ Loaded CD4: {len(df)} reactions, {df['significant'].sum()} significant")
        return df
    except FileNotFoundError:
        try:
            df = pd.read_csv('top_pathways_comprehensive.csv')
            st.sidebar.warning(f"⚠️ Using fallback CD4: {len(df)} reactions")
            return df
        except FileNotFoundError:
            st.sidebar.error("❌ CD4 data not found")
            return None

@st.cache_data
def load_cd8_data():
    """Load CD8 comprehensive pathway data"""
    try:
        df = pd.read_csv('cd8_pathways_comprehensive_RAW_PVALS.csv')
        st.sidebar.success(f"✅ Loaded CD8: {len(df)} reactions, {df['significant'].sum()} significant")
        return df
    except Exception as e:
        try:
            df = pd.read_csv('cd8_pathways_comprehensive_ttest.csv')
            st.sidebar.warning(f"⚠️ Using fallback CD8: {len(df)} reactions")
            return df
        except Exception as e2:
            st.sidebar.error(f"❌ CD8 data not found: {e2}")
            return None

@st.cache_data
def load_thymic_early_late():
    try:
        df = pd.read_csv('thymic_development_Early_vs_Late_Selection_comprehensive.csv')
        st.sidebar.success(f"✅ Early vs Late: {len(df)} reactions, {df['significant'].sum()} significant")
        return df
    except FileNotFoundError:
        st.sidebar.error("❌ Early vs Late data not found")
        return None

@st.cache_data
def load_thymic_late_mature():
    try:
        df = pd.read_csv('thymic_development_Late_vs_Mature_CD8SP_comprehensive.csv')
        st.sidebar.success(f"✅ Late vs Mature: {len(df)} reactions, {df['significant'].sum()} significant")
        return df
    except FileNotFoundError:
        st.sidebar.error("❌ Late vs Mature data not found")
        return None

@st.cache_data
def load_thymic_early_mature():
    try:
        df = pd.read_csv('thymic_development_Early_vs_Mature_CD8SP_comprehensive.csv')
        st.sidebar.success(f"✅ Early vs Mature: {len(df)} reactions, {df['significant'].sum()} significant")
        return df
    except FileNotFoundError:
        st.sidebar.error("❌ Early vs Mature data not found")
        return None

@st.cache_data
def load_thymic_three_way():
    try:
        # Try to load the comprehensive version with stage means first
        df = pd.read_csv('thymic_three_way_comprehensive_with_means.csv')
        st.sidebar.success(f"✅ Three-way ANOVA: {len(df)} reactions, {df['significant'].sum()} significant")
        return df
    except FileNotFoundError:
        try:
            # Fallback to original version
            df = pd.read_csv('thymic_three_way_anova_comprehensive.csv')
            st.sidebar.warning(f"⚠️ Using basic three-way data: {len(df)} reactions")
            return df
        except FileNotFoundError:
            st.sidebar.error("❌ Three-way comparison data not found")
            return None
    except Exception as e:
        st.sidebar.error(f"❌ Three-way data error: {e}")
        return None

def create_pathway_explorer(data, dataset_name, hi_label, lo_label):
    """Create the pathway explorer interface for a dataset"""
    
    if data is None:
        st.error(f"{dataset_name} data file not found. Please upload the data file.")
        return
    
    # Sidebar filters
    st.sidebar.header(f"🔍 {dataset_name} Search & Filters")
    
    # Search functionality
    st.sidebar.subheader("🔎 Search")
    search_term = st.sidebar.text_input(
        "Search pathways or genes:",
        placeholder="e.g., glycolysis, fatty acid, Pfkm, Ldha",
        help="Search pathway names and gene symbols",
        key=f"search_{dataset_name}"
    )
    
    st.sidebar.subheader("📊 Filters")
    
    # Direction filter by population/cell type
    direction_options = ['All'] + list(data['pathway_direction'].unique())
    selected_direction = st.sidebar.selectbox("Higher in Population", direction_options, key=f"direction_{dataset_name}")
    
    # Significance filter
    show_significant_only = st.sidebar.checkbox("Show only significant reactions (p < 0.05)", value=False, key=f"sig_{dataset_name}")
    
    # Gene filter
    show_with_genes_only = st.sidebar.checkbox("Show only reactions with known genes", value=False, key=f"genes_{dataset_name}")
    
    # Effect size threshold
    min_effect_size = st.sidebar.slider("Minimum |Effect Size|", 0.0, 5.0, 0.0, 0.1, key=f"effect_{dataset_name}")
    
    # Apply search filter
    filtered_data = data.copy()
    
    if search_term:
        search_term_lower = search_term.lower()
        
        # Search in pathway names and gene symbols automatically
        pathway_match = filtered_data['pathway'].str.lower().str.contains(search_term_lower, na=False)
        gene_match = filtered_data['genes'].str.lower().str.contains(search_term_lower, na=False)
        
        search_mask = pathway_match | gene_match
        filtered_data = filtered_data[search_mask]
        
        # Show search results summary
        if len(filtered_data) == 0:
            st.sidebar.warning(f"No results found for '{search_term}'")
        else:
            st.sidebar.success(f"Found {len(filtered_data)} reactions matching '{search_term}'")
    
    # Fix pathway_direction based on pathway median Cohen's d (for thymic data)
    if hi_label != lo_label:  # Two-group comparison
        filtered_data['pathway_direction_fixed'] = filtered_data['pathway_median_d'].apply(
            lambda x: hi_label if x > 0 else lo_label
        )
        direction_col = 'pathway_direction_fixed'
    else:  # Three-way comparison
        direction_col = 'pathway_direction'
    
    # Apply other filters
    if selected_direction != 'All':
        filtered_data = filtered_data[filtered_data[direction_col] == selected_direction]
    
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
        
        # Show search suggestions if no results
        if search_term and len(pathways) == 0:
            st.info("💡 **Try searching for:**")
            st.write("- **Pathways:** glycolysis, fatty acid, oxidative, transport")
            st.write("- **Genes:** Pfkm, Ldha, Cs, Idh1, Slc2a1")
        
        # Pathway selector with enhanced display
        selected_pathways = []
        for _, pathway_row in pathways.iterrows():
            pathway_name = pathway_row['pathway']
            
            # COLOR LOGIC: Red for positive Cohen's d, Blue for negative Cohen's d
            direction_emoji = "🔴" if pathway_row['median_d'] > 0 else "🔵"
            
            # Enhanced label with effect size and significance
            effect_size = pathway_row['median_d']
            pct_sig = pathway_row['pct_significant']
            
            # Truncate pathway name and add stats
            short_name = pathway_name[:35] + "..." if len(pathway_name) > 35 else pathway_name
            checkbox_label = f"{direction_emoji} {short_name} (d={effect_size:+.2f}, {pct_sig:.0f}% sig)"
            
            if st.checkbox(checkbox_label, key=f"pathway_{dataset_name}_{pathway_name}"):
                selected_pathways.append(pathway_name)
        
        # Show selected pathway summary
        if selected_pathways:
            st.subheader("📊 Selected Pathways")
            for pathway in selected_pathways:
                pathway_info = pathways[pathways['pathway'] == pathway].iloc[0]
                # COLOR LOGIC: Red for positive Cohen's d, Blue for negative Cohen's d
                direction_color = "🔴" if pathway_info['median_d'] > 0 else "🔵"
                st.write(f"{direction_color} **{pathway[:40]}...**")
                st.write(f"   Effect: {pathway_info['median_d']:+.2f} | Significant: {pathway_info['n_significant']}/{pathway_info['n_total']} ({pathway_info['pct_significant']:.1f}%)")
    
    with col2:
        st.header("🧪 Reactions")
        
        # Show search results summary in main area
        if search_term:
            st.info(f"🔍 Showing results for: **{search_term}**")
        
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
                format_func=lambda x: f"{x} (d={pathway_data[pathway_data['reaction_id']==x]['cohens_d'].iloc[0]:+.2f})",
                key=f"reaction_{dataset_name}"
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
    
def create_three_way_pathway_explorer(data, dataset_name):
    """Create specialized pathway explorer for three-way ANOVA data"""
    
    if data is None:
        st.error(f"{dataset_name} data file not found. Please upload the data file.")
        return
    
    # For three-way data, use Cohen's d and developmental trajectory
    st.info("🔬 **Three-way ANOVA analysis** - Showing developmental patterns across Early Selection, Late Selection, and Mature CD8SP stages")
    
    # Sidebar filters
    st.sidebar.header(f"🔍 {dataset_name} Search & Filters")
    
    # Search functionality
    st.sidebar.subheader("🔎 Search")
    search_term = st.sidebar.text_input(
        "Search pathways or genes:",
        placeholder="e.g., glycolysis, fatty acid, Pfkm, Ldha",
        help="Search pathway names and gene symbols",
        key=f"search_{dataset_name}"
    )
    
    st.sidebar.subheader("📊 Filters")
    
    # Trajectory filter
    if 'developmental_trajectory' in data.columns:
        trajectory_options = ['All'] + list(data['developmental_trajectory'].unique())
        selected_trajectory = st.sidebar.selectbox("Developmental Pattern", trajectory_options, key=f"trajectory_{dataset_name}")
    else:
        selected_trajectory = 'All'
        st.sidebar.warning("⚠️ Developmental trajectory data not available")
    
    # Significance filter
    show_significant_only = st.sidebar.checkbox("Show only significant reactions (p < 0.05)", value=False, key=f"sig_{dataset_name}")
    
    # Gene filter
    show_with_genes_only = st.sidebar.checkbox("Show only reactions with known genes", value=False, key=f"genes_{dataset_name}")
    
    # Effect size threshold
    min_effect_size = st.sidebar.slider("Minimum |Cohen's d|", 0.0, 5.0, 0.0, 0.1, key=f"effect_{dataset_name}")
    
    # Apply search filter
    filtered_data = data.copy()
    
    if search_term:
        search_term_lower = search_term.lower()
        
        # Search in pathway names and gene symbols automatically
        pathway_match = filtered_data['pathway'].str.lower().str.contains(search_term_lower, na=False)
        gene_match = filtered_data['genes'].str.lower().str.contains(search_term_lower, na=False)
        
        search_mask = pathway_match | gene_match
        filtered_data = filtered_data[search_mask]
        
        # Show search results summary
        if len(filtered_data) == 0:
            st.sidebar.warning(f"No results found for '{search_term}'")
        else:
            st.sidebar.success(f"Found {len(filtered_data)} reactions matching '{search_term}'")
    
    # Apply other filters
    if selected_trajectory != 'All' and 'developmental_trajectory' in filtered_data.columns:
        filtered_data = filtered_data[filtered_data['developmental_trajectory'] == selected_trajectory]
    
    if show_significant_only:
        filtered_data = filtered_data[filtered_data['significant'] == True]
    
    if show_with_genes_only:
        filtered_data = filtered_data[filtered_data['n_genes'] > 0]
    
    # Filter by effect size
    filtered_data = filtered_data[abs(filtered_data['cohens_d']) >= min_effect_size]
    
    # Main content area
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.header("📋 Pathways")
        
        # Get unique pathways from filtered data with stats
        agg_dict = {
            'cohens_d': 'median',
            'significant': 'sum',
            'reaction_id': 'count'
        }
        
        if 'developmental_trajectory' in filtered_data.columns:
            agg_dict['developmental_trajectory'] = lambda x: x.value_counts().index[0]  # Most common trajectory
            
        pathways = filtered_data.groupby('pathway').agg(agg_dict).reset_index()
        
        if 'developmental_trajectory' in pathways.columns:
            pathways.columns = ['pathway', 'median_d', 'n_significant', 'n_total', 'common_trajectory']
        else:
            pathways.columns = ['pathway', 'median_d', 'n_significant', 'n_total']
            pathways['common_trajectory'] = 'Unknown'
            
        pathways['pct_significant'] = (pathways['n_significant'] / pathways['n_total'] * 100).round(1)
        
        # Sort pathways by effect size (absolute value) then by % significant
        pathways['abs_median_d'] = abs(pathways['median_d'])
        pathways = pathways.sort_values(['abs_median_d', 'pct_significant'], ascending=[False, False])
        pathways = pathways.drop('abs_median_d', axis=1)
        
        st.write(f"**{len(pathways)} pathways match filters**")
        
        # Show search suggestions if no results
        if search_term and len(pathways) == 0:
            st.info("💡 **Try searching for:**")
            st.write("- **Pathways:** glycolysis, fatty acid, oxidative, transport")
            st.write("- **Genes:** Pfkm, Ldha, Cs, Idh1, Slc2a1")
        
        # Pathway selector with enhanced display for three-way data
        selected_pathways = []
        for _, pathway_row in pathways.iterrows():
            pathway_name = pathway_row['pathway']
            
            # COLOR LOGIC: Based on which stage has highest activity
            if 'highest_stage' in pathway_row and pd.notna(pathway_row['highest_stage']):
                # New comprehensive data with individual stage means
                highest_stage = pathway_row['highest_stage']
                if highest_stage == 'Mature_CD8SP':
                    stage_emoji = "🔴"  # Red for Mature highest
                elif highest_stage == 'Late_Selection':
                    stage_emoji = "🟡"  # Yellow for Late highest
                elif highest_stage == 'Early_Selection':
                    stage_emoji = "🔵"  # Blue for Early highest
                else:
                    stage_emoji = "⚪"  # White for unclear
            else:
                # Fallback for old data format
                effect_size = pathway_row['median_d']
                if effect_size > 0:
                    stage_emoji = "🔴"
                elif effect_size < 0:
                    stage_emoji = "🔵"
                else:
                    stage_emoji = "⚪"
            
            # Enhanced label with Cohen's d and trajectory
            effect_size = pathway_row['median_d']
            pct_sig = pathway_row['pct_significant']
            
            # Truncate pathway name and add stats
            short_name = pathway_name[:30] + "..." if len(pathway_name) > 30 else pathway_name
            checkbox_label = f"{stage_emoji} {short_name} (d={effect_size:+.2f}, {pct_sig:.0f}% sig)"
            
            if st.checkbox(checkbox_label, key=f"pathway_{dataset_name}_{pathway_name}"):
                selected_pathways.append(pathway_name)
        
        # Show selected pathway summary
        if selected_pathways:
            st.subheader("📊 Selected Pathways")
            for pathway in selected_pathways:
                pathway_info = pathways[pathways['pathway'] == pathway].iloc[0]
                
                # COLOR LOGIC: Based on which stage has highest activity
                if 'highest_stage' in pathway_info and pd.notna(pathway_info['highest_stage']):
                    # New comprehensive data with individual stage means
                    highest_stage = pathway_info['highest_stage']
                    if highest_stage == 'Mature_CD8SP':
                        stage_color = "🔴"  # Red for Mature highest
                    elif highest_stage == 'Late_Selection':
                        stage_color = "🟡"  # Yellow for Late highest
                    elif highest_stage == 'Early_Selection':
                        stage_color = "🔵"  # Blue for Early highest
                    else:
                        stage_color = "⚪"  # White for unclear
                    
                    # Show ranking with actual values if available
                    if 'ranking_text' in pathway_info and pd.notna(pathway_info['ranking_text']):
                        ranking_display = f"Ranking: {pathway_info['ranking_text']}"
                    else:
                        ranking_display = f"Highest: {highest_stage.replace('_', ' ')}"
                else:
                    # Fallback for old data format
                    effect_size = pathway_info['median_d']
                    if effect_size > 0:
                        stage_color = "🔴"
                    elif effect_size < 0:
                        stage_color = "🔵"
                    else:
                        stage_color = "⚪"
                    ranking_display = f"Pattern: {pathway_info['common_trajectory']} | Effect: {pathway_info['median_d']:+.2f}"
                
                st.write(f"{stage_color} **{pathway[:40]}...**")
                st.write(f"   {ranking_display} | Significant: {pathway_info['n_significant']}/{pathway_info['n_total']} ({pathway_info['pct_significant']:.1f}%)")
    
    with col2:
        st.header("🧪 Reactions")
        
        # Show search results summary in main area
        if search_term:
            st.info(f"🔍 Showing results for: **{search_term}**")
        
        if selected_pathways:
            # Filter data for selected pathways
            pathway_data = filtered_data[filtered_data['pathway'].isin(selected_pathways)]
            
            # Sort by effect size (absolute value)
            pathway_data = pathway_data.sort_values('cohens_d', key=abs, ascending=False)
            
            # Display reactions table
            st.write(f"**{len(pathway_data)} reactions in selected pathway(s)**")
            
            # Create display dataframe for three-way data
            display_columns = ['reaction_id', 'pathway', 'cohens_d', 'p_value', 'significant', 'genes', 'n_genes', 'reaction_name']
            if 'developmental_trajectory' in pathway_data.columns:
                display_columns.insert(5, 'developmental_trajectory')
                
            display_df = pathway_data[display_columns].copy()
            display_df['pathway_short'] = display_df['pathway'].str[:25] + "..."
            display_df['genes_short'] = display_df['genes'].str[:40] + "..."
            display_df['cohens_d'] = display_df['cohens_d'].round(2)
            display_df['p_value'] = display_df['p_value'].apply(lambda x: f"{x:.2e}")
            display_df['significant'] = display_df['significant'].map({True: "✓", False: "✗"})
            
            if 'developmental_trajectory' in display_df.columns:
                display_df['trajectory_short'] = display_df['developmental_trajectory'].str[:15] + "..."
                show_columns = ['reaction_id', 'pathway_short', 'cohens_d', 'p_value', 'significant', 'trajectory_short', 'n_genes', 'genes_short']
                column_config = {
                    "reaction_id": "Reaction ID",
                    "pathway_short": "Pathway",
                    "cohens_d": "Cohen's d",
                    "p_value": "p-value",
                    "significant": "Sig",
                    "trajectory_short": "Pattern",
                    "n_genes": "# Genes",
                    "genes_short": "Associated Genes"
                }
            else:
                show_columns = ['reaction_id', 'pathway_short', 'cohens_d', 'p_value', 'significant', 'n_genes', 'genes_short']
                column_config = {
                    "reaction_id": "Reaction ID",
                    "pathway_short": "Pathway",
                    "cohens_d": "Cohen's d",
                    "p_value": "p-value",
                    "significant": "Sig",
                    "n_genes": "# Genes",
                    "genes_short": "Associated Genes"
                }
            
            # Display table
            st.dataframe(
                display_df[show_columns],
                column_config=column_config,
                hide_index=True,
                use_container_width=True
            )
            
            # Reaction details for three-way data
            st.subheader("🔍 Reaction Details")
            selected_reaction = st.selectbox(
                "Select a reaction for details:",
                options=pathway_data['reaction_id'].tolist(),
                format_func=lambda x: f"{x} (d={pathway_data[pathway_data['reaction_id']==x]['cohens_d'].iloc[0]:+.2f})",
                key=f"reaction_{dataset_name}"
            )
            
            if selected_reaction:
                reaction_info = pathway_data[pathway_data['reaction_id'] == selected_reaction].iloc[0]
                
                col_a, col_b = st.columns(2)
                
                with col_a:
                    st.write(f"**Reaction:** {selected_reaction}")
                    st.write(f"**Name:** {reaction_info['reaction_name']}")
                    st.write(f"**Pathway:** {reaction_info['pathway']}")
                    st.write(f"**Cohen's d:** {reaction_info['cohens_d']:+.2f}")
                    st.write(f"**p-value:** {reaction_info['p_value']:.2e}")
                    st.write(f"**Significant:** {'Yes' if reaction_info['significant'] else 'No'}")
                    if 'developmental_trajectory' in reaction_info:
                        st.write(f"**Pattern:** {reaction_info['developmental_trajectory']}")
                
                with col_b:
                    # Show stage means if available (new comprehensive data)
                    if 'early_mean' in reaction_info and 'late_mean' in reaction_info and 'mature_mean' in reaction_info:
                        st.write(f"**Stage Activity Levels:**")
                        st.write(f"- Early Selection: {reaction_info['early_mean']:.2f}")
                        st.write(f"- Late Selection: {reaction_info['late_mean']:.2f}")
                        st.write(f"- Mature CD8SP: {reaction_info['mature_mean']:.2f}")
                        
                        if 'ranking_text' in reaction_info and pd.notna(reaction_info['ranking_text']):
                            st.write(f"**Ranking:** {reaction_info['ranking_text']}")
                    
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
    
    # Summary statistics for three-way data
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

# LOAD ALL DATA FIRST
cd4_data_loaded = load_cd4_data()
cd8_data_loaded = load_cd8_data()
thymic_early_late_loaded = load_thymic_early_late()
thymic_late_mature_loaded = load_thymic_late_mature()
thymic_early_mature_loaded = load_thymic_early_mature()
thymic_three_way_loaded = load_thymic_three_way()

# CREATE TABS
tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    " CD4+ T Cells", 
    " CD8+ T Cells", 
    " Early vs Late", 
    " Late vs Mature", 
    " Early vs Mature",
    " Three-Way ANOVA"
])

with tab1:
    st.header("CD4+ T Cell Metabolic Activity")
    st.markdown("""
    **Explore differential metabolic flux between CD5 hi and CD5 lo CD4+ T cells**
    - Data from Compass metabolic flux analysis
    - **ALL pathways** ordered by effect size and significance
    - Gene associations from Mouse-GEM metabolic model
    - 🔴 **Red means positive Cohen's d (CD5 hi > CD5 lo)**, 🔵 **Blue means negative Cohen's d (CD5 lo > CD5 hi)**
    - Cohen's D measures 'how different' two groups are (0 = no difference)
    - Cohen's D of +2 = CD5 hi cells have higher metabolic flux, -2 = CD5 lo cells have higher flux
    """)
    create_pathway_explorer(cd4_data_loaded, "CD4", "CD5 hi", "CD5 lo")

with tab2:
    st.header("CD8+ T Cell Metabolic Activity")
    st.markdown("""
    **Explore differential metabolic flux between CD5 hi and CD5 lo CD8+ T cells**
    - Data from Compass metabolic flux analysis
    - **ALL pathways** ordered by effect size and significance
    - Gene associations from Mouse-GEM metabolic model
    - 🔴 **Red means positive Cohen's d (CD8 CD5 hi > CD8 CD5 lo)**, 🔵 **Blue means negative Cohen's d (CD8 CD5 lo > CD8 CD5 hi)**
    - Cohen's D measures 'how different' two groups are (0 = no difference)
    - Cohen's D of +2 = CD8 CD5 hi cells have higher metabolic flux, -2 = CD8 CD5 lo cells have higher flux
    """)
    create_pathway_explorer(cd8_data_loaded, "CD8", "CD8_CD5_hi", "CD8_CD5_lo")

with tab3:
    st.header("Thymic Development: Early vs Late Selection")
    st.markdown("""
    **Explore metabolic changes during T cell positive selection**
    - Comparison: Early Selection vs Late Selection stages
    - Data from Compass metabolic flux analysis of thymic development
    - Gene associations from Mouse-GEM metabolic model
    - 🔴 **Red means positive Cohen's d (Early Selection > Late Selection)**, 🔵 **Blue means negative Cohen's d (Late Selection > Early Selection)**
    - Cohen's D measures 'how different' the developmental stages are (0 = no difference)
    - Cohen's D of +2 = Early Selection higher flux, -2 = Late Selection higher flux
    """)
    create_pathway_explorer(thymic_early_late_loaded, "Early_vs_Late", "Early_Selection", "Late_Selection")

with tab4:
    st.header("Thymic Development: Late Selection vs Mature CD8SP")
    st.markdown("""
    **Explore metabolic changes during T cell maturation**
    - Comparison: Late Selection vs Mature CD8SP stages
    - Data from Compass metabolic flux analysis of thymic development
    - Gene associations from Mouse-GEM metabolic model
    - 🔴 **Red means positive Cohen's d (Late Selection > Mature CD8SP)**, 🔵 **Blue means negative Cohen's d (Mature CD8SP > Late Selection)**
    - Cohen's D measures 'how different' the developmental stages are (0 = no difference)
    - Cohen's D of +2 = Mature CD8SP higher flux, -2 = Late Selection higher flux
    """)
    create_pathway_explorer(thymic_late_mature_loaded, "Late_vs_Mature", "Late_Selection", "Mature_CD8SP")

with tab5:
    st.header("Thymic Development: Early Selection vs Mature CD8SP")
    st.markdown("""
    **Explore metabolic changes across T cell development**
    - Comparison: Early Selection vs Mature CD8SP stages
    - Data from Compass metabolic flux analysis of thymic development
    - Gene associations from Mouse-GEM metabolic model
    - 🔴 **Red means positive Cohen's d (Early Selection > Mature CD8SP)**, 🔵 **Blue means negative Cohen's d (Mature CD8SP > Early Selection)**
    - Cohen's D measures 'how different' the developmental stages are (0 = no difference)
    - Cohen's D of +2 = Mature CD8SP higher flux, -2 = Early Selection higher flux
    """)
    create_pathway_explorer(thymic_early_mature_loaded, "Early_vs_Mature", "Early_Selection", "Mature_CD8SP")

with tab6:
    st.header("Thymic Development: Three-Way ANOVA Comparison")
    st.markdown("""
    **Explore reactions that differ significantly across all three developmental stages**
    - ANOVA analysis: Early Selection, Late Selection, and Mature CD8SP
    - Data from Compass metabolic flux analysis of thymic development
    - Gene associations from Mouse-GEM metabolic model
    - 🔴 **Red = Highest in Mature CD8SP**, 🟡 **Yellow = Highest in Late Selection**, 🔵 **Blue = Highest in Early Selection**, ⚪ **White = Unclear**
    - Shows which developmental stage has the highest metabolic activity for each pathway
    - Ranking text shows exact activity levels: "Late Selection (2.5) > Early Selection (1.8) > Mature CD8SP (1.2)"
    """)
    if thymic_three_way_loaded is not None:
        create_three_way_pathway_explorer(thymic_three_way_loaded, "Three_Way_ANOVA")
    else:
        st.error("Three-way comparison data not available yet")
        st.info("💡 Try clicking the '🔄 Clear Cache' button in the sidebar and refresh the page if the file exists.")
