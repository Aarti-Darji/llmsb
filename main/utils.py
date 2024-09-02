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
