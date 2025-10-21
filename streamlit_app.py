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
        df = pd.read_csv('thymic_three_way_anova_comprehensive.csv')
        st.sidebar.success(f"✅ Three-way ANOVA: {len(df)} reactions, {df['significant'].sum()} significant")
        return df
    except FileNotFoundError:
        st.sidebar.error("❌ Three-way comparison data not found")
        return None

# Replace the create_pathway_explorer function in your streamlit_app.py with this corrected version:

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
    
    # Direction filter
    direction_options = ['All'] + list(data['pathway_direction'].unique())
    selected_direction = st.sidebar.selectbox("Pathway Direction", direction_options, key=f"direction_{dataset_name}")
    
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
    
    # Apply other filters
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
        
        # Show search suggestions if no results
        if search_term and len(pathways) == 0:
            st.info("💡 **Try searching for:**")
            st.write("- **Pathways:** glycolysis, fatty acid, oxidative, transport")
            st.write("- **Genes:** Pfkm, Ldha, Cs, Idh1, Slc2a1")
        
        # Pathway selector with enhanced display
        selected_pathways = []
        for _, pathway_row in pathways.iterrows():
            pathway_name = pathway_row['pathway']
            
            # FIXED COLOR LOGIC: Use exact string comparison, not substring matching
            direction_emoji = "🔴" if pathway_row['direction'] == hi_label else "🔵"
            
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
                # FIXED COLOR LOGIC: Use exact string comparison
                direction_color = "🔴" if pathway_info['direction'] == hi_label else "🔵"
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
    
    # Trajectory filter (instead of direction)
    if 'developmental_trajectory' in data.columns:
        trajectory_options = ['All'] + list(data['developmental_trajectory'].unique())
        selected_trajectory = st.sidebar.selectbox("Developmental Pattern", trajectory_options, key=f"trajectory_{dataset_name}")
    else:
        selected_trajectory = 'All'
        st.sidebar.warning("⚠️ Developmental trajectory data not available")
    
    # Stage with highest activity filter - for three-way data, we need to determine this from the means
    if 'early_mean' in data.columns and 'late_mean' in data.columns and 'mature_mean' in data.columns:
        # Calculate highest stage for each reaction
        data_copy = data.copy()
        def get_highest_stage(row):
            means = {
                'Early_Selection': row['early_mean'],
                'Late_Selection': row['late_mean'], 
                'Mature_CD8SP': row['mature_mean']
            }
            return max(means, key=means.get)
        
        data_copy['highest_stage'] = data_copy.apply(get_highest_stage, axis=1)
        stage_options = ['All'] + list(data_copy['highest_stage'].unique())
        selected_stage = st.sidebar.selectbox("Highest Activity Stage", stage_options, key=f"stage_{dataset_name}")
    else:
        stage_options = ['All'] + list(data['pathway_direction'].unique())
        selected_stage = st.sidebar.selectbox("Highest Activity Stage", stage_options, key=f"stage_{dataset_name}")
    
    # Significance filter
    show_significant_only = st.sidebar.checkbox("Show only significant reactions (p < 0.05)", value=False, key=f"sig_{dataset_name}")
    
    # Gene filter
    show_with_genes_only = st.sidebar.checkbox("Show only reactions with known genes", value=False, key=f"genes_{dataset_name}")
    
    # F-statistic threshold (instead of Cohen's d)
    min_f_stat = st.sidebar.slider("Minimum F-statistic", 0.0, 50.0, 0.0, 1.0, key=f"f_stat_{dataset_name}")
    
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
    
    if selected_stage != 'All':
        if 'early_mean' in filtered_data.columns and 'late_mean' in filtered_data.columns and 'mature_mean' in filtered_data.columns:
            # For three-way data, filter by calculated highest stage
            def get_highest_stage(row):
                means = {
                    'Early_Selection': row['early_mean'],
                    'Late_Selection': row['late_mean'], 
                    'Mature_CD8SP': row['mature_mean']
                }
                return max(means, key=means.get)
            
            filtered_data['highest_stage'] = filtered_data.apply(get_highest_stage, axis=1)
            filtered_data = filtered_data[filtered_data['highest_stage'] == selected_stage]
        else:
            filtered_data = filtered_data[filtered_data['pathway_direction'] == selected_stage]
    
    if show_significant_only:
        filtered_data = filtered_data[filtered_data['significant'] == True]
    
    if show_with_genes_only:
        filtered_data = filtered_data[filtered_data['n_genes'] > 0]
    
    # Filter by F-statistic
    filtered_data = filtered_data[filtered_data['f_statistic'] >= min_f_stat]
    
    # Main content area
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.header("📋 Pathways")
        
        # Get unique pathways from filtered data with stats
        if 'early_mean' in filtered_data.columns and 'late_mean' in filtered_data.columns and 'mature_mean' in filtered_data.columns:
            # For three-way data, calculate highest stage for each pathway
            if 'highest_stage' not in filtered_data.columns:
                def get_highest_stage(row):
                    means = {
                        'Early_Selection': row['early_mean'],
                        'Late_Selection': row['late_mean'], 
                        'Mature_CD8SP': row['mature_mean']
                    }
                    return max(means, key=means.get)
                filtered_data['highest_stage'] = filtered_data.apply(get_highest_stage, axis=1)
            
            agg_dict = {
                'highest_stage': lambda x: x.value_counts().index[0],  # Most common highest stage
                'f_statistic': 'median',
                'significant': 'sum',
                'reaction_id': 'count'
            }
            
            if 'developmental_trajectory' in filtered_data.columns:
                agg_dict['developmental_trajectory'] = lambda x: x.value_counts().index[0]  # Most common trajectory
                
            pathways = filtered_data.groupby('pathway').agg(agg_dict).reset_index()
            
            if 'developmental_trajectory' in pathways.columns:
                pathways.columns = ['pathway', 'highest_stage', 'median_f', 'n_significant', 'n_total', 'common_trajectory']
            else:
                pathways.columns = ['pathway', 'highest_stage', 'median_f', 'n_significant', 'n_total']
                pathways['common_trajectory'] = 'Unknown'
        else:
            # Regular two-group comparison
            pathways = filtered_data.groupby('pathway').agg({
                'pathway_direction': 'first',
                'developmental_trajectory': lambda x: x.value_counts().index[0] if 'developmental_trajectory' in filtered_data.columns else 'N/A',
                'f_statistic': 'median',
                'significant': 'sum',
                'reaction_id': 'count'
            }).reset_index()
            
            pathways.columns = ['pathway', 'highest_stage', 'common_trajectory', 'median_f', 'n_significant', 'n_total']
        pathways['pct_significant'] = (pathways['n_significant'] / pathways['n_total'] * 100).round(1)
        
        # Sort pathways by F-statistic then by % significant
        pathways = pathways.sort_values(['median_f', 'pct_significant'], ascending=[False, False])
        
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
            
            # Color by highest stage
            if pathway_row['highest_stage'] == 'Mature_CD8SP':
                stage_emoji = "🔴"
            elif pathway_row['highest_stage'] == 'Early_Selection':
                stage_emoji = "🔵"
            else:
                stage_emoji = "🟡"  # Late Selection
            
            # Enhanced label with F-statistic and trajectory
            f_stat = pathway_row['median_f']
            pct_sig = pathway_row['pct_significant']
            trajectory = pathway_row['common_trajectory']
            
            # Truncate pathway name and add stats
            short_name = pathway_name[:30] + "..." if len(pathway_name) > 30 else pathway_name
            checkbox_label = f"{stage_emoji} {short_name} (F={f_stat:.1f}, {pct_sig:.0f}% sig)"
            
            if st.checkbox(checkbox_label, key=f"pathway_{dataset_name}_{pathway_name}"):
                selected_pathways.append(pathway_name)
        
        # Show selected pathway summary
        if selected_pathways:
            st.subheader("📊 Selected Pathways")
            for pathway in selected_pathways:
                pathway_info = pathways[pathways['pathway'] == pathway].iloc[0]
                if pathway_info['highest_stage'] == 'Mature_CD8SP':
                    stage_color = "🔴"
                elif pathway_info['highest_stage'] == 'Early_Selection':
                    stage_color = "🔵"
                else:
                    stage_color = "🟡"
                st.write(f"{stage_color} **{pathway[:40]}...**")
                st.write(f"   Highest: {pathway_info['highest_stage']} | F-stat: {pathway_info['median_f']:.1f} | Significant: {pathway_info['n_significant']}/{pathway_info['n_total']} ({pathway_info['pct_significant']:.1f}%)")
    
    with col2:
        st.header("🧪 Reactions")
        
        # Show search results summary in main area
        if search_term:
            st.info(f"🔍 Showing results for: **{search_term}**")
        
        if selected_pathways:
            # Filter data for selected pathways
            pathway_data = filtered_data[filtered_data['pathway'].isin(selected_pathways)]
            
            # Sort by F-statistic
            pathway_data = pathway_data.sort_values('f_statistic', ascending=False)
            
            # Display reactions table
            st.write(f"**{len(pathway_data)} reactions in selected pathway(s)**")
            
            # Create display dataframe for three-way data
            display_df = pathway_data[['reaction_id', 'pathway', 'f_statistic', 'p_value', 'significant', 'developmental_trajectory', 'pathway_direction', 'genes', 'n_genes', 'reaction_name']].copy()
            display_df['pathway_short'] = display_df['pathway'].str[:25] + "..."
            display_df['genes_short'] = display_df['genes'].str[:40] + "..."
            display_df['f_statistic'] = display_df['f_statistic'].round(2)
            display_df['p_value'] = display_df['p_value'].apply(lambda x: f"{x:.2e}")
            display_df['significant'] = display_df['significant'].map({True: "✓", False: "✗"})
            display_df['trajectory_short'] = display_df['developmental_trajectory'].str[:15] + "..."
            
            # Select columns to show
            show_columns = ['reaction_id', 'pathway_short', 'f_statistic', 'p_value', 'significant', 'pathway_direction', 'trajectory_short', 'n_genes', 'genes_short']
            
            # Display table
            st.dataframe(
                display_df[show_columns],
                column_config={
                    "reaction_id": "Reaction ID",
                    "pathway_short": "Pathway",
                    "f_statistic": "F-statistic",
                    "p_value": "p-value",
                    "significant": "Sig",
                    "pathway_direction": "Highest Stage",
                    "trajectory_short": "Pattern",
                    "n_genes": "# Genes",
                    "genes_short": "Associated Genes"
                },
                hide_index=True,
                use_container_width=True
            )
            
            # Reaction details for three-way data
            st.subheader("🔍 Reaction Details")
            selected_reaction = st.selectbox(
                "Select a reaction for details:",
                options=pathway_data['reaction_id'].tolist(),
                format_func=lambda x: f"{x} (F={pathway_data[pathway_data['reaction_id']==x]['f_statistic'].iloc[0]:.1f})",
                key=f"reaction_{dataset_name}"
            )
            
            if selected_reaction:
                reaction_info = pathway_data[pathway_data['reaction_id'] == selected_reaction].iloc[0]
                
                col_a, col_b = st.columns(2)
                
                with col_a:
                    st.write(f"**Reaction:** {selected_reaction}")
                    st.write(f"**Name:** {reaction_info['reaction_name']}")
                    st.write(f"**Pathway:** {reaction_info['pathway']}")
                    st.write(f"**F-statistic:** {reaction_info['f_statistic']:.2f}")
                    st.write(f"**p-value:** {reaction_info['p_value']:.2e}")
                    st.write(f"**Significant:** {'Yes' if reaction_info['significant'] else 'No'}")
                    st.write(f"**Highest Stage:** {reaction_info['pathway_direction']}")
                    st.write(f"**Pattern:** {reaction_info['developmental_trajectory']}")
                
                with col_b:
                    st.write(f"**Stage Means:**")
                    st.write(f"- Early Selection: {reaction_info['early_mean']:.2f}")
                    st.write(f"- Late Selection: {reaction_info['late_mean']:.2f}")
                    st.write(f"- Mature CD8SP: {reaction_info['mature_mean']:.2f}")
                    
                    st.write(f"**Associated Genes ({reaction_info['n_genes']}):**")
                    if reaction_info['genes'] and reaction_info['genes'] != '':
                        genes_list = reaction_info['genes'].split('; ')
                        for gene in genes_list:
                            st.write(f"- {gene}")
                    else:
                        st.write("No genes associated")
        
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
    "Three-Way ANOVA"
])

with tab1:
    st.header("CD4+ T Cell Metabolic Activity")
    st.markdown("""
    **Explore differential metabolic flux between CD5 hi and CD5 lo CD4+ T cells**
    - Data from Compass metabolic flux analysis
    - **ALL pathways** ordered by effect size and significance
    - Gene associations from Mouse-GEM metabolic model
    - 🔴 **Red means higher in CD5 hi**, 🔵 **Blue means higher in CD5 lo**
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
    - 🔴 **Red means higher in CD8 CD5 hi**, 🔵 **Blue means higher in CD8 CD5 lo**
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
    - 🔴 **Red means higher in Early Selection**, 🔵 **Blue means higher in Late Selection**
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
    - 🔴 **Red means higher in Mature CD8SP**, 🔵 **Blue means higher in Late Selection**
    - Cohen's D measures 'how different' the developmental stages are (0 = no difference)
    - Cohen's D of +2 = Mature CD8SP higher flux, -2 = Late Selection higher flux
    """)
    create_pathway_explorer(thymic_late_mature_loaded, "Late_vs_Mature", "Mature_CD8SP", "Late_Selection")

with tab5:
    st.header("Thymic Development: Early Selection vs Mature CD8SP")
    st.markdown("""
    **Explore metabolic changes across T cell development**
    - Comparison: Early Selection vs Mature CD8SP stages
    - Data from Compass metabolic flux analysis of thymic development
    - Gene associations from Mouse-GEM metabolic model
    - 🔴 **Red means higher in Mature CD8SP**, 🔵 **Blue means higher in Early Selection**
    - Cohen's D measures 'how different' the developmental stages are (0 = no difference)
    - Cohen's D of +2 = Mature CD8SP higher flux, -2 = Early Selection higher flux
    """)
    create_pathway_explorer(thymic_early_mature_loaded, "Early_vs_Mature", "Mature_CD8SP", "Early_Selection")

with tab6:
    st.header("Thymic Development: Three-Way ANOVA Comparison")
    st.markdown("""
    **Explore reactions that differ significantly across all three developmental stages**
    - ANOVA analysis: Early Selection, Late Selection, and Mature CD8SP
    - Data from Compass metabolic flux analysis of thymic development
    - Gene associations from Mouse-GEM metabolic model
    - 🔴 **Red = Highest in Mature CD8SP**, 🟡 **Yellow = Highest in Late Selection**, 🔵 **Blue = Highest in Early Selection**
    - F-statistic measures overall difference across all three stages
    - Developmental patterns: Increasing, Decreasing, Peak, Dip, Complex
    """)
    if thymic_three_way_loaded is not None:
        st.info("🔬 **Three-way ANOVA analysis** showing reactions that differ significantly across Early Selection, Late Selection, and Mature CD8SP stages")
        create_three_way_pathway_explorer(thymic_three_way_loaded, "Three_Way_ANOVA")
    else:
        st.error("Three-way comparison data not available yet")
