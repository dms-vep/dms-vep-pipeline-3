<script>
import * as vega from 'vega';
import * as vegaLite from 'vega-lite';
import vegaEmbed from 'vega-embed';
import { parseVegaSpecFromHTML } from './parseVegaSpec';

export default {
    name: 'VegaChart',

    props: {
        specUrl: {
            type: String,
            required: true,
        },
        showShadow: {
            type: Boolean,
            default: true,
        }
    },

    data() {
        return {
            isExpanded: false,
        };
    },

    async mounted() {
        await this.loadChart(this.specUrl);
    },

    methods: {
        // Parse and load the chart from the URL
        async loadChart(url) {
            try {
                let spec;
                // Determine if the URL is local or remote
                if (!url.startsWith('http://') && !url.startsWith('https://')) {
                    const htmlContent = await this.fetchLocalFile(url);
                    spec = await parseVegaSpecFromHTML(htmlContent);
                } else {

                    // Get the path to check the file format
                    const urlObj = new URL(url);
                    const path = urlObj.pathname;
                    // Get the response from the URL
                    const response = await fetch(url);

                    // Check if URL ends with '.html'
                    if (path.endsWith('.html')) {
                        const htmlContent = await response.text();
                        spec = await parseVegaSpecFromHTML(htmlContent);
                    }
                    // Check if URL ends with '.json'
                    else if (path.endsWith('.json')) {
                        spec = await response.json();
                    }
                    else {
                        console.error('Unsupported file format');
                        return;
                    }
                }
                this.renderChart(spec);
            } catch (error) {
                console.error('Error loading Vega spec:', error);
            }
        },
        // Fetch the text from local files
        async fetchLocalFile(filePath) {
            const response = await fetch(filePath);
            if (!response.ok) throw new Error('Failed to fetch local file');
            console.log('response.text()', response);
            return await response.text();
        },
        // Render the chart using VegaEmbed
        renderChart(spec) {
            vegaEmbed(this.$refs.vegaContainer, spec, {
                renderer: 'canvas',
                vega,
                vegaLite,
                actions: false,
            }).then((result) => {
                console.log('Chart rendered successfully');
            }).catch(console.error);
        },
        // Toggle the expanded view
        toggleExpand() {
            this.isExpanded = !this.isExpanded;
        },
    },
};
</script>

<template>
    <div class="vega-chart-container" :class="{ 'expanded': isExpanded, 'no-box-shadow': !showShadow }">
        <div ref="vegaContainer">
        </div>
        <button @click="toggleExpand" @keydown.esc="isExpanded = false" class="expand-btn">
            <i v-if="isExpanded" class="bi bi-arrows-angle-contract"></i>
            <i v-else class="bi bi-arrows-angle-expand"></i>
        </button>
    </div>
</template>

<style>
.vega-chart-container {
    width: 100%;
    height: auto;
    overflow: auto;
    padding: 20px;
    position: relative;
    background-color: white;
    line-height: 0%;
}

.vega-chart-container:not(.no-box-shadow) {
    box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.5);
}

.expanded {
    width: 100vw;
    height: 100vh;
    position: fixed !important;
    top: 0;
    left: 0;
    z-index: 100;
    padding: 1%;
}

.expand-btn {
    position: absolute;
    top: 10px;
    right: 10px;
    outline: none !important;
    width: 30px;
    height: 30px;
    border-radius: 5px;
    border: 1px solid #aaa;
    display: flex;
    justify-content: center;
    align-items: center;
}

.bi {
    font-size: 20px;
}

/* Form container */
.vega-bindings {
    display: flex;
    flex-direction: column;
    gap: 10px;
    padding: 10px;
}

.vega-bind-name {
    flex-basis: 20%;
    text-align: left;
    font-weight: bold;
    margin-right: 5px;
}

.vega-bindings select {
    border: 1px solid #ccc;
    border-radius: 4px;
    background-color: #fff;
    text-align: center;
}

.vega-bindings input[type="select"] option {
    text-align: center;
    font-size: 14px;
}

.vega-bindings input[type="radio"] {
    margin: 5px 5px 0px 10px;
}
</style>