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
def load_strain_early():
    try:
        df = pd.read_csv('strain_comparison_early_selection_comprehensive.csv')
        st.sidebar.success(f"✅ Early Strain: {len(df)} reactions, {df['significant'].sum()} significant")
        return df
    except FileNotFoundError:
        st.sidebar.error("❌ Early strain comparison data not found")
        return None

@st.cache_data
def load_strain_late():
    try:
        df = pd.read_csv('strain_comparison_late_selection_comprehensive.csv')
        st.sidebar.success(f"✅ Late Strain: {len(df)} reactions, {df['significant'].sum()} significant")
        return df
    except FileNotFoundError:
        st.sidebar.error("❌ Late strain comparison data not found")
        return None

@st.cache_data
def load_strain_mature():
    try:
        df = pd.read_csv('strain_comparison_mature_cd8sp_comprehensive.csv')
        st.sidebar.success(f"✅ Mature Strain: {len(df)} reactions, {df['significant'].sum()} significant")
        return df
    except FileNotFoundError:
        st.sidebar.error("❌ Mature strain comparison data not found")
        return None

def create_strain_pathway_explorer(data, dataset_name, stage_name):
    """Create the pathway explorer interface for strain comparisons (F5 vs OT1 vs TG6)"""
    
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
    
    # Strain filter (instead of direction)
    strain_options = ['All'] + list(data['highest_strain'].unique())
    selected_strain = st.sidebar.selectbox("Highest Activity Strain", strain_options, key=f"strain_{dataset_name}")
    
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
    
    # Apply strain filter - need to use pathway-level aggregation
    if selected_strain != 'All':
        # Get pathways where the selected strain is highest
        pathway_strain_map = filtered_data.groupby('pathway')['pathway_highest_strain'].first()
        pathways_for_strain = pathway_strain_map[pathway_strain_map == selected_strain].index
        filtered_data = filtered_data[filtered_data['pathway'].isin(pathways_for_strain)]
    
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
            'pathway_highest_strain': 'first',
            'pathway_median_d': 'first',
            'significant': 'sum',
            'reaction_id': 'count'
        }).reset_index()
        
        pathways.columns = ['pathway', 'highest_strain', 'median_d', 'n_significant', 'n_total']
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
            
            # Color based on highest strain (by CD5 expression level)
            strain_colors = {"OT1": "🔴", "F5": "🟡", "TG6": "🔵"}
            strain_emoji = strain_colors.get(pathway_row['highest_strain'], "⚪")
            
            # Enhanced label with effect size and significance
            effect_size = pathway_row['median_d']
            pct_sig = pathway_row['pct_significant']
            
            # Truncate pathway name and add stats
            short_name = pathway_name[:35] + "..." if len(pathway_name) > 35 else pathway_name
            checkbox_label = f"{strain_emoji} {short_name} (d={effect_size:+.2f}, {pct_sig:.0f}% sig)"
            
            if st.checkbox(checkbox_label, key=f"pathway_{dataset_name}_{pathway_name}"):
                selected_pathways.append(pathway_name)
        
        # Show selected pathway summary
        if selected_pathways:
            st.subheader("📊 Selected Pathways")
            for pathway in selected_pathways:
                pathway_info = pathways[pathways['pathway'] == pathway].iloc[0]
                strain_colors = {"OT1": "🔴", "F5": "🟡", "TG6": "🔵"}
                strain_color = strain_colors.get(pathway_info['highest_strain'], "⚪")
                st.write(f"{strain_color} **{pathway[:40]}...**")
                st.write(f"   Highest: {pathway_info['highest_strain']} | Effect: {pathway_info['median_d']:+.2f} | Significant: {pathway_info['n_significant']}/{pathway_info['n_total']} ({pathway_info['pct_significant']:.1f}%)")
    
    with col2:
        st.header("🧪 Reactions")
        
        if selected_pathways:
            # Filter data for selected pathways
            pathway_data = filtered_data[filtered_data['pathway'].isin(selected_pathways)]
            
            if len(pathway_data) > 0:
                # Show reaction details
                st.write(f"**{len(pathway_data)} reactions in selected pathways**")
                
                # Show detailed table
                display_cols = ['reaction_id', 'pathway', 'f_statistic', 'p_value', 'cohens_d', 
                               'ranking_text', 'highest_strain', 'reaction_name', 'genes']
                
                if len(pathway_data) <= 1000:
                    st.dataframe(
                        pathway_data[display_cols].round(4),
                        use_container_width=True,
                        height=400
                    )
                else:
                    st.warning(f"Showing first 1000 of {len(pathway_data)} reactions")
                    st.dataframe(
                        pathway_data[display_cols].head(1000).round(4),
                        use_container_width=True,
                        height=400
                    )
            else:
                st.info("No reactions found for selected pathways with current filters.")
        else:
            st.info("👈 Select pathways from the sidebar to view detailed reaction data")
    
    # Summary statistics at bottom
    st.subheader("📊 Dataset Summary")
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Total Pathways", len(data['pathway'].unique()))
    
    with col2:
        st.metric("Total Reactions", len(data))
    
    with col3:
        st.metric("Significant Reactions", len(data[data['significant'] == True]))
    
    with col4:
        st.metric("Reactions with Genes", len(data[data['n_genes'] > 0]))

def create_pathway_explorer(data, dataset_name, hi_label, lo_label):
    """Create the pathway explorer interface for two-group comparisons"""
    
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
            
            # Color based on Cohen's d sign (red for positive, blue for negative)
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
                direction_color = "🔴" if pathway_info['median_d'] > 0 else "🔵"
                st.write(f"{direction_color} **{pathway[:40]}...**")
                st.write(f"   Effect: {pathway_info['median_d']:+.2f} | Significant: {pathway_info['n_significant']}/{pathway_info['n_total']} ({pathway_info['pct_significant']:.1f}%)")
    
    with col2:
        st.header("🧪 Reactions")
        
        if selected_pathways:
            # Filter data for selected pathways
            pathway_data = filtered_data[filtered_data['pathway'].isin(selected_pathways)]
            
            if len(pathway_data) > 0:
                # Show reaction details
                st.write(f"**{len(pathway_data)} reactions in selected pathways**")
                
                # Show detailed table
                display_cols = ['reaction_id', 'pathway', 'cohens_d', 'p_value', 'significant', 
                               'reaction_name', 'genes']
                
                if len(pathway_data) <= 1000:
                    st.dataframe(
                        pathway_data[display_cols].round(4),
                        use_container_width=True,
                        height=400
                    )
                else:
                    st.warning(f"Showing first 1000 of {len(pathway_data)} reactions")
                    st.dataframe(
                        pathway_data[display_cols].head(1000).round(4),
                        use_container_width=True,
                        height=400
                    )
        else:
                st.info("No reactions found for selected pathways with current filters.")
    
    # Summary statistics at bottom
    st.subheader("📊 Dataset Summary")
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
strain_early_loaded = load_strain_early()
strain_late_loaded = load_strain_late()
strain_mature_loaded = load_strain_mature()

# CREATE TABS
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "CD4+ T Cells", 
    "CD8+ T Cells", 
    "Early Strains",
    "Late Strains", 
    "Mature Strains"
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
    st.header("Strain Comparison: Early Selection Stage")
    st.markdown("""
    **Explore metabolic differences between F5, OT1, and TG6 strains during early T cell selection**
    - ANOVA analysis: F5 vs OT1 vs TG6 at Early Selection stage
    - Data from Compass metabolic flux analysis of thymic development
    - Gene associations from Mouse-GEM metabolic model
    - 🔴 **Red = Highest in OT1 (High CD5)**, 🟡 **Yellow = Highest in F5 (Intermediate CD5)**, 🔵 **Blue = Highest in TG6 (Low CD5)**
    - Shows which strain has the highest metabolic activity for each pathway
    - Ranking text shows exact activity levels: "TG6 (8761.57) > F5 (8733.29) > OT1 (8650.87)"
    """)
    create_strain_pathway_explorer(strain_early_loaded, "Early_Strains", "Early Selection")

with tab4:
    st.header("Strain Comparison: Late Selection Stage")
    st.markdown("""
    **Explore metabolic differences between F5, OT1, and TG6 strains during late T cell selection**
    - ANOVA analysis: F5 vs OT1 vs TG6 at Late Selection stage
    - Data from Compass metabolic flux analysis of thymic development
    - Gene associations from Mouse-GEM metabolic model
    - 🔴 **Red = Highest in OT1 (High CD5)**, 🟡 **Yellow = Highest in F5 (Intermediate CD5)**, 🔵 **Blue = Highest in TG6 (Low CD5)**
    - Shows which strain has the highest metabolic activity for each pathway
    - Ranking text shows exact activity levels for each reaction
    """)
    create_strain_pathway_explorer(strain_late_loaded, "Late_Strains", "Late Selection")

with tab5:
    st.header("Strain Comparison: Mature CD8SP Stage")
    st.markdown("""
    **Explore metabolic differences between F5, OT1, and TG6 strains in mature CD8+ T cells**
    - ANOVA analysis: F5 vs OT1 vs TG6 at Mature CD8SP stage
    - Data from Compass metabolic flux analysis of thymic development
    - Gene associations from Mouse-GEM metabolic model
    - 🔴 **Red = Highest in OT1 (High CD5)**, 🟡 **Yellow = Highest in F5 (Intermediate CD5)**, 🔵 **Blue = Highest in TG6 (Low CD5)**
    - Shows which strain has the highest metabolic activity for each pathway
    - Ranking text shows exact activity levels for each reaction
    """)
    create_strain_pathway_explorer(strain_mature_loaded, "Mature_Strains", "Mature CD8SP")

# Add Cohen's d explanation at the bottom
st.markdown("---")
st.subheader("📖 Understanding Cohen's d")
st.markdown("""
**Cohen's d** is a standardized measure of effect size that quantifies the difference between two groups:

- **Cohen's d = 0**: No difference between groups
- **Cohen's d = ±0.2**: Small effect (subtle difference)
- **Cohen's d = ±0.5**: Medium effect (moderate difference)  
- **Cohen's d = ±0.8**: Large effect (substantial difference)
- **Cohen's d = ±2.0**: Very large effect (major difference)

**Interpretation for metabolic data:**
- **Positive Cohen's d**: First group has higher metabolic flux
- **Negative Cohen's d**: Second group has higher metabolic flux
- **Larger absolute value**: Greater difference in metabolic activity

For strain comparisons, we show which strain has the highest activity and provide ranking text with exact values.
""")
