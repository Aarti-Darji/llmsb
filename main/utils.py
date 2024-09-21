from langchain_core.tools import tool
import duckdb
import numpy as np

@tool
def find_high_interaction_regions(channel_ids, top_n):
    """
    Finds regions with high concentrations of the specified channels.

    Args:
        channel_ids (list): A list of channel IDs to search for high concentration regions.
        top_n (int): The number of top regions to return based on concentration levels.

    Returns:
        list: A list of tuples containing the (x, y) coordinates of the top regions.
    """
    con = duckdb.connect("/workspaces/llmsb/data/Dataset1-LSP13626-melanoma-in-situ-256_bricks (1).db")
    ids = [i + 1 for i in channel_ids]
    results = {channel_id: np.zeros((22, 43)) for channel_id in ids}

    for channel_id in ids:
        query = f"""
        SELECT bricks.chunk_key,
               SPLIT_PART(bricks.chunk_key, '.', 5) AS x,
               SPLIT_PART(bricks.chunk_key, '.', 4) AS y,
               bricks.mean
        FROM bricks
        WHERE bricks.channel_id = {channel_id}
          AND SPLIT_PART(bricks.chunk_key, '.', 3) = '100'
          AND bricks.mean > (
            SELECT AVG(mean)
            FROM bricks AS brick_b
            WHERE bricks.channel_id = brick_b.channel_id
              AND brick_b.mean > 0
            GROUP BY brick_b.channel_id
        );
        """

        res = con.execute(query)
        queryresult = res.fetchall()

        for entry in queryresult:
            x = int(entry[1])
            y = int(entry[2])
            mean = entry[3]
            results[channel_id][y, x] = mean

        if np.max(results[channel_id]) != 0:
            results[channel_id] = results[channel_id] / np.max(results[channel_id])

    combined_result = np.zeros((22, 43))
    for channel_id in ids:
        combined_result += results[channel_id]

    flat_combined_result = combined_result.flatten()
    top_indices = np.argsort(flat_combined_result)[-top_n:][::-1]
    top_coords = np.unravel_index(top_indices, combined_result.shape)
    top_coords = list(zip(top_coords[0], top_coords[1]))

    final_coords = []
    for (y, x) in top_coords:
        y = int(y) + 0.5
        x = int(x) + 0.5
        x = x * 256
        y = y * 256
        final_coords.append((x, y))

    con.close()
    return final_coords

@tool
def find_low_interaction_regions(channel_ids, top_n):
    """
    Finds regions with low concentrations of the specified channels.

    Args:
        channel_ids (list): A list of channel IDs to search for low concentration regions.
        top_n (int): The number of top regions to return based on low concentration levels.

    Returns:
        list: A list of tuples containing the (x, y) coordinates of the top low-concentration regions.
    """
    con = duckdb.connect("/workspaces/llmsb/data/Dataset1-LSP13626-melanoma-in-situ-256_bricks (1).db")
    ids = [i + 1 for i in channel_ids]
    results = {channel_id: np.zeros((22, 43)) for channel_id in ids}

    for channel_id in ids:
        query = f"""
        SELECT bricks.chunk_key,
               SPLIT_PART(bricks.chunk_key, '.', 5) AS x,
               SPLIT_PART(bricks.chunk_key, '.', 4) AS y,
               bricks.mean
        FROM bricks
        WHERE bricks.channel_id = {channel_id}
          AND SPLIT_PART(bricks.chunk_key, '.', 3) = '100'
          AND bricks.mean < (
            SELECT AVG(mean)
            FROM bricks AS brick_b
            WHERE bricks.channel_id = brick_b.channel_id
              AND brick_b.mean > 0
            GROUP BY brick_b.channel_id
        );
        """

        res = con.execute(query)
        queryresult = res.fetchall()

        for entry in queryresult:
            x = int(entry[1])
            y = int(entry[2])
            mean = entry[3]
            results[channel_id][y, x] = mean

        if np.max(results[channel_id]) != 0:
            results[channel_id] = results[channel_id] / np.max(results[channel_id])

    combined_result = np.zeros((22, 43))
    for channel_id in ids:
        combined_result += results[channel_id]

    flat_combined_result = combined_result.flatten()
    top_indices = np.argsort(flat_combined_result)[:top_n]  # Get indices for the lowest values
    top_coords = np.unravel_index(top_indices, combined_result.shape)
    top_coords = list(zip(top_coords[0], top_coords[1]))

    final_coords = []
    for (y, x) in top_coords:
        y = int(y) + 0.5
        x = int(x) + 0.5
        x = x * 256
        y = y * 256
        final_coords.append((x, y))

    con.close()
    return final_coords

@tool
def find_high_low_interaction_regions(high_marker_channel_ids, low_marker_channel_ids, top_n):
    """
    Finds regions with high concentrations of specified high markers and low concentrations of specified low markers.

    Args:
        high_marker_channel_ids (list): A list of channel IDs to search for high concentration regions.
        low_marker_channel_ids (list): A list of channel IDs to search for low concentration regions.
        top_n (int): The number of top regions to return based on high and low concentration levels.

    Returns:
        list: A list of tuples containing the (x, y) coordinates of the top regions with high and low concentrations.
    """
    con = duckdb.connect("/workspaces/llmsb/data/Dataset1-LSP13626-melanoma-in-situ-256_bricks (1).db")

    high_results = {}
    low_results = {}

    for channel_id in high_marker_channel_ids:
        query = f"""
        SELECT SPLIT_PART(chunk_key, '.', 5) AS x,
               SPLIT_PART(chunk_key, '.', 4) AS y,
               mean
        FROM bricks
        WHERE channel_id = {channel_id}
          AND SPLIT_PART(chunk_key, '.', 3) = '100'
          AND mean > (
            SELECT AVG(mean)
            FROM bricks
            WHERE channel_id = {channel_id}
              AND mean > 0
          )
        """
        res = con.execute(query)
        queryresult = res.fetchall()

        results = np.zeros((22, 43))
        for entry in queryresult:
            x = int(entry[0])
            y = int(entry[1])
            mean = entry[2]
            results[y, x] = mean

        if np.max(results) != 0:
            results = results / np.max(results)

        high_results[channel_id] = results

    for channel_id in low_marker_channel_ids:
        query = f"""
        SELECT SPLIT_PART(chunk_key, '.', 5) AS x,
               SPLIT_PART(chunk_key, '.', 4) AS y,
               mean
        FROM bricks
        WHERE channel_id = {channel_id}
          AND SPLIT_PART(chunk_key, '.', 3) = '100'
          AND mean < (
            SELECT AVG(mean)
            FROM bricks
            WHERE channel_id = {channel_id}
              AND mean > 0
          )
        """
        res = con.execute(query)
        queryresult = res.fetchall()

        results = np.zeros((22, 43))
        for entry in queryresult:
            x = int(entry[0])
            y = int(entry[1])
            mean = entry[2]
            results[y, x] = mean

        if np.max(results) != 0:
            results = results / np.max(results)

        low_results[channel_id] = results

    combined_high_result = np.min([high_results[channel_id] for channel_id in high_marker_channel_ids], axis=0)

    combined_low_result = np.max([low_results[channel_id] for channel_id in low_marker_channel_ids], axis=0)

    final_result = combined_high_result - combined_low_result

    flat_combined_result = final_result.flatten()
    top_indices = np.argsort(flat_combined_result)[-top_n:][::-1]
    top_coords = np.unravel_index(top_indices, final_result.shape)
    top_coords = list(zip(top_coords[0], top_coords[1]))

    final_coords = []
    for (y, x) in top_coords:
        y = int(y) + 0.5
        x = int(x) + 0.5
        x = x * 256
        y = y * 256
        final_coords.append((x, y))

    con.close()
    return final_coords

@tool
def find_channel_means_at_coords(coords: list) -> dict:
    """
    Finds the mean values of channels at the specified coordinates.

    Args:
        coords (list): A list of (x, y) tuples representing the coordinates to find the channel means.

    Returns:
        dict: A dictionary where keys are channel IDs and values are tuples containing
              (combined region mean, average mean across bricks, maximum mean, minimum mean).
    """
    con = duckdb.connect("/workspaces/llmsb/data/Dataset1-LSP13626-melanoma-in-situ-256_bricks (1).db")

    channel_means = {}

    for (x_coord, y_coord) in coords:
        x_chunk = x_coord / 256
        y_chunk = y_coord / 256
        y_chunk = int(y_chunk)
        x_chunk = int(x_chunk)
        query = f"""
        SELECT bricks.channel_id, bricks.mean
        FROM bricks
        WHERE SPLIT_PART(bricks.chunk_key, '.', 5) = '{x_chunk}'
          AND SPLIT_PART(bricks.chunk_key, '.', 4) = '{y_chunk}'
          AND SPLIT_PART(bricks.chunk_key, '.', 3) = '100';
        """
        res = con.execute(query)
        queryresult = res.fetchall()
        for entry in queryresult:
            channel_id = entry[0]
            brick_mean = entry[1]
            if brick_mean is None or brick_mean <= 0:
                continue  

            if channel_id not in channel_means:
                channel_means[channel_id] = []

            channel_means[channel_id].append(brick_mean)

    final_channel_means = {}

    for channel_id, brick_means in channel_means.items():
        if not brick_means:
            continue  

        combined_region_mean = sum(brick_means) / len(brick_means)

        avg_query = f"""
        SELECT AVG(mean)
        FROM bricks
        WHERE channel_id = {channel_id}
          AND mean > 0;
        """
        avg_res = con.execute(avg_query)
        avg_mean = avg_res.fetchone()[0]

        max_query = f"""
        SELECT MAX(mean)
        FROM bricks
        WHERE channel_id = {channel_id};
        """
        max_res = con.execute(max_query)
        max_mean = max_res.fetchone()[0]

        min_query = f"""
        SELECT MIN(mean)
        FROM bricks
        WHERE channel_id = {channel_id};
        """
        min_res = con.execute(min_query)
        min_mean = min_res.fetchone()[0]

        adjusted_channel_id = channel_id - 1
        final_channel_means[adjusted_channel_id] = (combined_region_mean, avg_mean, max_mean, min_mean)

    con.close()
    return final_channel_means

def invoke_agent(message):
    payload = {"input": message}
    responseagent = agent_executor.invoke(payload)
    output_text = responseagent['output']
    descr = output_text.split("```python")[0].strip()
    return descr, output_text

def extract_brick_ids(output_text):
    brick_ids_match = re.search(r'brick_ids\s*=\s*\[(.*?)\]', output_text)
    if brick_ids_match:
        brick_ids_str = brick_ids_match.group(1)
        brick_ids = [brick_id.strip().strip('"').strip("'") for brick_id in brick_ids_str.split(',')]
        return brick_ids
    return []

def update_config_with_brick_ids(old_config, brick_ids_to_highlight):
    selection_name = "Highlighted Bricks"
    coordination_space = old_config['coordinationSpace']

    if 'additionalObsSets' not in coordination_space:
        coordination_space['additionalObsSets'] = {}

    if 'A' not in coordination_space['additionalObsSets'] or coordination_space['additionalObsSets']['A'] is None:
        coordination_space['additionalObsSets']['A'] = {'version': '0.1.3', 'datatype': 'obs', 'tree': []}

    tree_structure = coordination_space['additionalObsSets']['A'].get('tree', [])
    my_selections = next((item for item in tree_structure if item['name'] == 'My Selections'), None)

    if my_selections is None:
        my_selections = {'name': 'My Selections', 'children': []}
        tree_structure.append(my_selections)

    highlighted_selection = next((child for child in my_selections['children'] if child['name'] == selection_name), None)
    if highlighted_selection is not None:
        highlighted_selection['set'] = [[brick_id, None] for brick_id in brick_ids_to_highlight]
    else:
        my_selections['children'].append({
            'name': selection_name,
            'set': [[brick_id, None] for brick_id in brick_ids_to_highlight]
        })

    coordination_space['additionalObsSets']['A']['tree'] = tree_structure
    coordination_space['obsSetSelection']['A'] = [['My Selections', selection_name]]

    if 'obsSetColor' not in coordination_space:
        coordination_space['obsSetColor'] = {'A': []}

    color_assignment = next((item for item in coordination_space['obsSetColor']['A'] if item['path'] == ['My Selections', selection_name]), None)
    if color_assignment is not None:
        color_assignment['color'] = [0, 0, 255]
    else:
        coordination_space['obsSetColor']['A'].append({
            'path': ['My Selections', selection_name],
            'color': [0, 0, 255]
        })

    coordination_space.setdefault('obsColorEncoding', {})
    coordination_space['obsColorEncoding']['init_bv_obsSegmentations_0'] = 'cellSetSelection'
    coordination_space['obsColorEncoding']['A'] = 'cellSetSelection'

    old_config["uid"] = f"with_spatial_target_{uuid.uuid4().hex}"

    return old_config

def parse_and_merge_config(old_config, output_text):
    code_block_match = re.search(r'```python\n(.*?)\n```', output_text, re.DOTALL)
    if not code_block_match:
        print("No code block found in agent's response.")
        return old_config  # Return the old configuration if no match is found

    code_block = code_block_match.group(1).strip()
    code_block = re.sub(r'CL\(\s*\[', '[', code_block)
    code_block = re.sub(r'\]\s*\)', ']', code_block)
    code_block = code_block.replace('True', 'true').replace('False', 'false')

    if not (code_block.startswith('{') and code_block.endswith('}')):
        code_block = '{' + code_block + '}'

    try:
        config_dict = json.loads(code_block)
    except json.JSONDecodeError as e:
        print(f"Error parsing code block as JSON: {e}")
        return old_config

    coordination_space = old_config.get('coordinationSpace', {})

    if 'imageLayer' in config_dict and isinstance(config_dict['imageLayer'], list):
        image_channels = config_dict['imageLayer']

        for idx, channel in enumerate(image_channels):
            if idx >= 6:
                break  # Only process up to 6 channels (init_bv_image_0 to init_bv_image_5)

            if 'spatialTargetC' in channel:
                channel_id = f'init_bv_image_{idx}'
                coordination_space.setdefault('spatialTargetC', {})[channel_id] = channel['spatialTargetC']
        for i in range(6):
            channel_id = f'init_bv_image_{i}'
            if i < len(image_channels):
                # If the channel is in the returned list, make it visible
                coordination_space.setdefault('spatialChannelVisible', {})[channel_id] = True
            else:
                # If the channel is not in the returned list, make it invisible
                coordination_space.setdefault('spatialChannelVisible', {})[channel_id] = False

    old_config['coordinationSpace'] = coordination_space

    old_config['uid'] = f"with_spatial_target_{uuid.uuid4().hex}"

    return old_config

def apply_config_to_widget(vw, new_config):
    vw.config = new_config

def extract_selected_bricks(config):
    selection_ids = []

    if (config.get('coordinationSpace') and
        config['coordinationSpace'].get('obsSetSelection') and
        config['coordinationSpace']['obsSetSelection'].get('A') and
        config['coordinationSpace']['obsSetSelection']['A']):

        selection_path = config['coordinationSpace']['obsSetSelection']['A'][0]

        additional_obs_sets = config['coordinationSpace'].get('additionalObsSets', {}).get('A')

        if additional_obs_sets and 'tree' in additional_obs_sets:
            additional_obs_sets_tree = additional_obs_sets['tree']

            selection_details = None
            for item in additional_obs_sets_tree:
                if item['name'] == selection_path[0]:
                    for child in item.get('children', []):
                        if child['name'] == selection_path[1]:
                            selection_details = child['set']

            if selection_details:
                # No conversion to int, since selection IDs are now in 'y.x' format
                selection_ids = [item[0] for item in selection_details]  # Leave as string
                print(selection_details)
    return selection_ids
