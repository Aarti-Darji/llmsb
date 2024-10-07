@tool
def find_high_interaction_regions(channel_ids, top_n):
    """Locate high-concentrated regions for the given channel_ids and return the top coordinates in y.x format."""

    import duckdb
    import numpy as np

    con = duckdb.connect('/content/drive/MyDrive/Dataset1-LSP13626-melanoma-in-situ-256_bricks.db')
    ids = [i + 1 for i in channel_ids]
    results = {channel_id: np.zeros((22, 43)) for channel_id in ids}

    # Step 1: Execute the SQL query and gather results for each requested channel
    for channel_id in ids:
        query = f"""
        SELECT bricks.id,
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

        # Step 2: Collect mean values for each coordinate
        for entry in queryresult:
            x = int(entry[1])
            y = int(entry[2])
            mean = entry[3]

            # Check array bounds
            if 0 <= y < 22 and 0 <= x < 43:
                results[channel_id][y, x] = mean
            else:
                print(f"Warning: x or y out of bounds. x={x}, y={y}")

        # Normalize the result for this channel
        if np.max(results[channel_id]) != 0:
            results[channel_id] = results[channel_id] / np.max(results[channel_id])

    # Step 3: Combine results across channels to find top regions
    combined_result = np.zeros((22, 43))
    for channel_id in ids:
        combined_result += results[channel_id]

    flat_combined_result = combined_result.flatten()
    top_indices = np.argsort(flat_combined_result)[-top_n:][::-1]
    top_coords = np.unravel_index(top_indices, combined_result.shape)
    top_coords = list(zip(top_coords[0], top_coords[1]))

    # Return the final coordinates in y.x format
    final_coords = [f"{y}_{x}" for y, x in top_coords]

    con.close()
    return final_coords

@tool
def find_low_interaction_regions(channel_ids, top_n):
    """Locate low-concentrated regions for the given channel_ids and return the top coordinates in y.x format."""

    con = duckdb.connect('/content/drive/MyDrive/Dataset1-LSP13626-melanoma-in-situ-256_bricks.db')
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

    # Flatten and sort to find regions with the lowest mean values
    flat_combined_result = combined_result.flatten()
    top_indices = np.argsort(flat_combined_result)[:top_n]  # Get indices for the lowest values
    top_coords = np.unravel_index(top_indices, combined_result.shape)
    top_coords = list(zip(top_coords[0], top_coords[1]))

    # Return the final coordinates in y.x format
    final_coords = [f"{y}_{x}" for y, x in top_coords]

    con.close()
    return final_coords

@tool
def find_high_low_interaction_regions(high_marker_channel_ids, low_marker_channel_ids, top_n):
    """Locate regions with high concentrations of several markers and low concentrations of other markers."""
    con = duckdb.connect('/content/drive/MyDrive/Dataset1-LSP13626-melanoma-in-situ-256_bricks.db')

    # Initialize result arrays for high and low markers
    high_results = {}
    low_results = {}

    # Query for high markers
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

    # Query for low markers
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

    # Combine results for high markers
    combined_high_result = np.min([high_results[channel_id] for channel_id in high_marker_channel_ids], axis=0)

    # Combine results for low markers
    combined_low_result = np.max([low_results[channel_id] for channel_id in low_marker_channel_ids], axis=0)

    # Calculate final result
    final_result = combined_high_result - combined_low_result

    # Get top N coordinates
    flat_combined_result = final_result.flatten()
    top_indices = np.argsort(flat_combined_result)[-top_n:][::-1]
    top_coords = np.unravel_index(top_indices, final_result.shape)
    top_coords = list(zip(top_coords[0], top_coords[1]))

    # Return the final coordinates in y.x format
    final_coords = [f"{y}_{x}" for y, x in top_coords]

    con.close()
    return final_coords

@tool
def find_channel_means_at_coords(coords: list, channel_ids: list = None) -> dict:
    """
    Locate the bricks with the given coordinates in 'y.x' format, and optionally filter by channel_ids.
    Add 1 to channel IDs provided.
    Return the brick mean, overall tissue mean, max mean, and min mean for each channel.
    If there are 6 or fewer channels and 10 or fewer coords, return results by bricks with per-channel statistics.
    Otherwise, return results by channels with combined statistics.

    Args:
    coords (list of str): A list of coordinates in 'y.x' format as strings.
    channel_ids (list of int, optional): A list of channel IDs to filter by. If None, return statistics for all channels.

    Returns:
    dict: A dictionary with either:
      - bricks as keys (if fewer than 6 channels and 10 coords), and per-channel statistics (mean, min, max)
      - channels as keys, and combined statistics (mean, min, max) for the region.
    """
    con = duckdb.connect('/content/drive/MyDrive/Dataset1-LSP13626-melanoma-in-situ-256_bricks.db')

    # Prepare a list of conditions for the query using the coords
    coord_conditions = []
    for coord in coords:
        y_chunk, x_chunk = coord.split('_')
        coord_conditions.append(f"(SPLIT_PART(bricks.chunk_key, '.', 5) = '{x_chunk}' AND SPLIT_PART(bricks.chunk_key, '.', 4) = '{y_chunk}')")

    # Join all the conditions with OR to form a single query
    coord_conditions_str = " OR ".join(coord_conditions)

    # Add 1 to each channel ID if provided
    if channel_ids:
        channel_ids = [channel_id + 1 for channel_id in channel_ids]  # Increment each channel_id by 1
        channel_ids_str = ", ".join(map(str, channel_ids))  # Convert channel IDs to a string format for SQL
        channel_condition = f"AND bricks.channel_id IN ({channel_ids_str})"
    else:
        channel_condition = ""  # No filtering by channel_id

    # Query to get brick mean values by coordinates and channels
    query = f"""
    SELECT SPLIT_PART(bricks.chunk_key, '.', 4) || '.' || SPLIT_PART(bricks.chunk_key, '.', 5) AS brick_coord,
           bricks.channel_id, bricks.mean
    FROM bricks
    WHERE ({coord_conditions_str})
      AND SPLIT_PART(bricks.chunk_key, '.', 3) = '100'
      {channel_condition};
    """

    # Execute the query and fetch results
    res = con.execute(query)
    queryresult = res.fetchall()

    # Check if we are in the "brick-wise" mode or not
    if len(coords) <= 10 and (channel_ids is None or len(channel_ids) <= 6):
        # Brick-wise structure
        brick_means = {}

        for entry in queryresult:
            brick_coord = entry[0]  # e.g., '14.11'
            channel_id = entry[1]
            brick_mean = entry[2]

            # If the brick isn't in the dictionary yet, initialize it
            if brick_coord not in brick_means:
                brick_means[brick_coord] = {}

            # Store the channel data under each brick
            if channel_id not in brick_means[brick_coord]:
                brick_means[brick_coord][channel_id] = {"mean": brick_mean}

        # Now get the min and max values for each channel across the tissue
        for brick in brick_means:
            for channel_id in brick_means[brick]:
                # Get min and max for the channel
                min_query = f"""
                SELECT MIN(mean)
                FROM bricks
                WHERE channel_id = {channel_id};
                """
                min_res = con.execute(min_query)
                min_mean = min_res.fetchone()[0]

                max_query = f"""
                SELECT MAX(mean)
                FROM bricks
                WHERE channel_id = {channel_id};
                """
                max_res = con.execute(max_query)
                max_mean = max_res.fetchone()[0]

                # Update the channel data with min and max
                brick_means[brick][channel_id].update({
                    "min": min_mean,
                    "max": max_mean
                })

        con.close()
        return brick_means
    else:
        # Channel-wise structure (original structure)
        channel_means = {}

        for entry in queryresult:
            channel_id = entry[1]
            brick_mean = entry[2]

            if channel_id not in channel_means:
                channel_means[channel_id] = []

            channel_means[channel_id].append(brick_mean)

        final_channel_means = {}
        for channel_id, brick_means in channel_means.items():
            if not brick_means:
                continue  # Skip channels that have no valid means

            combined_region_mean = sum(brick_means) / len(brick_means)

            # Query to get the average mean for the channel
            avg_query = f"""
            SELECT AVG(mean)
            FROM bricks
            WHERE channel_id = {channel_id}
              AND mean > 0;
            """
            avg_res = con.execute(avg_query)
            avg_mean = avg_res.fetchone()[0]

            # Query to get the max mean for the channel
            max_query = f"""
            SELECT MAX(mean)
            FROM bricks
            WHERE channel_id = {channel_id};
            """
            max_res = con.execute(max_query)
            max_mean = max_res.fetchone()[0]

            # Query to get the min mean for the channel
            min_query = f"""
            SELECT MIN(mean)
            FROM bricks
            WHERE channel_id = {channel_id};
            """
            min_res = con.execute(min_query)
            min_mean = min_res.fetchone()[0]

            adjusted_channel_id = channel_id - 1  # Adjust back for the final result
            final_channel_means[adjusted_channel_id] = (combined_region_mean, avg_mean, max_mean, min_mean)

        con.close()
        return final_channel_means

from langchain.schema import AIMessage, HumanMessage

def invoke_agent(message):
    inputs = {"messages": [("user", message)]}
    config = {"configurable": {"thread_id": "10"}}  # Use a fixed thread_id if needed
    output_text = ""
    try:
        stream = graph.stream(inputs, config=config, stream_mode="values")
        last_processed_id = None

        for s in stream:
            print(f"Stream output: {s}")
            message_data = s.get("messages", [])
            if not message_data:
                continue
            last_message = message_data[-1]
            if last_message.id == last_processed_id:
                continue

            last_processed_id = last_message.id

            if isinstance(last_message, AIMessage):
                output_text += last_message.content + "\n"
                print(f"AI Response: {last_message.content}")

    except Exception as e:
        print(f"Error occurred during agent invocation: {e}")
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
        # print("No code block found in agent's response.")
        return old_config

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
    """
    Extract the latest selected bricks (coordinates) and the associated channels/markers from the configuration.

    Args:
    config (dict): The vw.config dictionary.

    Returns:
    tuple: A tuple containing two lists:
           - A list of selection IDs (coordinates in 'y.x' format).
           - A list of channel IDs (markers) corresponding to those coordinates.
    """
    selection_info = []
    channels = []

    # Extract selection details (coordinates)
    if (config.get('coordinationSpace') and
        config['coordinationSpace'].get('obsSetSelection') and
        config['coordinationSpace']['obsSetSelection'].get('A') and
        config['coordinationSpace']['obsSetSelection']['A']):

        # Get the most recent selection (last entry in the list)
        selection_path = config['coordinationSpace']['obsSetSelection']['A'][-1]

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
                # Extract selection coordinates (in 'y.x' format)
                selection_info = [item[0] for item in selection_details]

        channelsFromSettings = vw.config['coordinationSpace']['metaCoordinationScopesBy']['init_bv_image_0']['imageLayer']['imageChannel']['init_bv_image_0']
        channels = []
        for key in channelsFromSettings:
            if(vw.config['coordinationSpace']['spatialChannelVisible'][key]):
                channels.append(vw.config['coordinationSpace']['spatialTargetC'][key]+1)


    return selection_info, channels


def send_message(message, buffers):
    global vw
    if vw is None:
        print("Widget not initialized.")
        return

    old_config = copy.deepcopy(vw.config)
    selected_bricks, channels = extract_selected_bricks(old_config)

    if selected_bricks:
        message += f" bricks: {selected_bricks} and following markers are currently highlighted: {channels}"
        # print(message)
    descr, output_text = invoke_agent(message)
    brick_ids_to_highlight = extract_brick_ids(output_text)
    new_config = parse_and_merge_config(old_config, output_text)

    if brick_ids_to_highlight:
        new_config = update_config_with_brick_ids(new_config, brick_ids_to_highlight)

    apply_config_to_widget(vw, new_config)

    return {**new_config, "text": descr}, []

class ChatPlugin(VitesscePlugin):
    plugin_esm = PLUGIN_ESM
    commands = {
        "chat_send": send_message,
    }